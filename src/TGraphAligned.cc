#include "TGraphAligned.h" 

#include <stdlib.h>
#include "TMath.h"
#include <stdlib.h>
static int err = 0; // error checking for posix_memalign
#include "TH1.h"
#include "TList.h"


ClassImp(TGraphAligned); 

#define ALIGNMENT TGRAPH_ALIGNED_ALIGNMENT

TGraphAligned::TGraphAligned()
  :TGraph()
{
  SetEditable(0); 
  //don't think we need to do anything in this case. 
}

TGraphAligned::TGraphAligned(Int_t n)
{
  SetTitle("Graph"); 
  SetName("Graph"); 
  fNpoints =n; 
  if (!CtorAllocate()) return; 
  FillZero(0, fNpoints); 
  SetEditable(0); 
}

TGraphAligned& TGraphAligned::operator=(const TGraphAligned &gr)
{
   // Equal operator for this graph

   if (this != &gr) {
      TNamed::operator=(gr);
      TAttLine::operator=(gr);
      TAttFill::operator=(gr);
      TAttMarker::operator=(gr);

      fNpoints = gr.fNpoints;
      fMaxSize = gr.fMaxSize;

      // delete list of functions and their contents before copying it
      if (fFunctions) {
         // delete previous lists of functions
         if (!fFunctions->IsEmpty()) {
            fFunctions->SetBit(kInvalidObject);
            // use TList::Remove to take into account the case the same object is
            // added multiple times in the list
            TObject *obj;
            while ((obj  = fFunctions->First())) {
               while (fFunctions->Remove(obj)) { }
               delete obj;
            }
         }
         delete fFunctions;
      }

      if (gr.fFunctions) fFunctions = (TList*)gr.fFunctions->Clone();
      else fFunctions = new TList;

      if (fHistogram) delete fHistogram;
      if (gr.fHistogram) fHistogram = new TH1F(*(gr.fHistogram));
      else fHistogram = 0;

      fMinimum = gr.fMinimum;
      fMaximum = gr.fMaximum;
      if (fX) free(fX); 
      if (fY) free(fY); 
      if (!fMaxSize) {
         fX = fY = 0;
         return *this;
      }
      else
      {
        err = posix_memalign((void**) &fX, ALIGNMENT, fNpoints * sizeof(Double_t));
        if(err != 0){
          fX = NULL;
        }
        err = posix_memalign((void**) &fY, ALIGNMENT, fNpoints * sizeof(Double_t));
        if(err != 0){
          fX = NULL;
        }	
      }

      Int_t n = gr.GetN() * sizeof(Double_t);
      if (n > 0) {
         memcpy(fX, gr.fX, n);
         memcpy(fY, gr.fY, n);
      }
   }
   SetEditable(0); 
   return *this;
}

TGraphAligned::TGraphAligned(Int_t n, const Double_t * x, const Double_t * y)
{

  SetTitle("Graph"); 
  SetName("Graph"); 
  if (!x || !y) {
      fNpoints = 0;
   } else {
      fNpoints = n;
   }
   if (!CtorAllocate()) return;
   n = fNpoints * sizeof(Double_t);
   memcpy(fX, x, n);
   memcpy(fY, y, n);

  SetEditable(0); 
}




// I don't think this needs to do everything other than allocate fX and fY since the TGraph constructor should take care of that
Bool_t TGraphAligned::CtorAllocate() 
{
  if (fNpoints <= 0)
  {
       fNpoints = 0;
       fMaxSize   = 0;
       fX         = 0;
       fY         = 0;
       return kFALSE;
  } 
  else 
  {
    fMaxSize   = fNpoints;
    err = posix_memalign((void**) &fX, ALIGNMENT, fNpoints * sizeof(Double_t));
    if(err != 0){
      fX = NULL;
    }

    err = posix_memalign((void**) &fY, ALIGNMENT, fNpoints * sizeof(Double_t));
    if(err != 0){
      fY = NULL;
     }
  }


  if (fX == NULL || fY == NULL) 
  {
    fprintf(stderr, "Could not allocate aligned memory in TGraphAligned::CtorAllocate()\n"); 
    return kFALSE; 
  }

  return kTRUE;
}

//same as TGraph except use free instead of delete[] 
TGraphAligned::~TGraphAligned()
{
   // Graph default destructor.

   free(fX);
   free(fY);
   fX = 0; 
   fY = 0; 
}


//same as TGraph except use free instead of delete[] 
void TGraphAligned::CopyAndRelease(Double_t **newarrays, Int_t ibegin, Int_t iend,
                            Int_t obegin)
{
   CopyPoints(newarrays, ibegin, iend, obegin);
   if (newarrays) {
      free(fX);
      fX = newarrays[0];
      free(fY);
      fY = newarrays[1];
      delete[] newarrays;
   }
}

TGraphAligned::TGraphAligned(const TGraph &gr)
   : TGraph(gr) 
{
   // Copy constructor for this graph

  // Ok , we'll let the TGraph copy constructor do all the work but then 
  // remake the pointers. If anybody can think of a better way to do this... let me know
   if (fMaxSize)
   {

     Double_t * new_fX = NULL;
     err = posix_memalign((void**) &new_fX, ALIGNMENT, fMaxSize * sizeof(Double_t));
     if(err != 0){
       new_fX = NULL;
     }
     Double_t * new_fY = NULL;
     err = posix_memalign((void**) &new_fY, ALIGNMENT, fMaxSize * sizeof(Double_t));
     if(err != 0){
       new_fY = NULL;
     }
     
     memcpy(new_fX, fX, gr.GetN() * sizeof(Double_t));
     delete [] fX; 
     fX = new_fX; 


     memcpy(new_fY, fY, gr.GetN() * sizeof(Double_t));
     delete [] fY; 
     fY = new_fY; 
   }
}


//same as TGraph except use memalign instead of new[] 
Double_t** TGraphAligned::AllocateAlignedArrays(Int_t Narrays, Int_t arraySize)
{
  if (arraySize < 0)
  {
    arraySize = 0;
  }

  Double_t **newarrays = new Double_t*[Narrays];
  if (!arraySize)
  {
    for (Int_t i = 0; i < Narrays; ++i)
    {
       newarrays[i] = 0;
    }
  } 
  else
  {
    for (Int_t i = 0; i < Narrays; ++i)
    {
      err = posix_memalign((void**) &newarrays[i], ALIGNMENT, arraySize * sizeof(Double_t));
      if(err != 0)
      {
        newarrays[i] = NULL;
      }
    }
  }

  fMaxSize = arraySize;
  return newarrays;
}


void TGraphAligned::undBize() 
{
  int n = GetN(); 

  for (int i = 0; i < n; i++)
  {
    fY[i] = TMath::Power(10,fY[i]/10); 
  }
}


Double_t TGraphAligned::getSumV2(Int_t istart, Int_t iend) const
{
  aligned_double_v v = GetY(); 
  __builtin_prefetch(v);  // This is really an academic exercise at this point
  int N = GetN();  
  int real_start = istart >= 0 ? TMath::Min(istart, N-1) : TMath::Max(0, N + istart)  ; 
  int real_end = iend >= 0 ? TMath::Min(iend, N-1) : TMath::Max(0, N + iend)  ; 

  double sum2 = 0; 


  for (int i = real_start; i <= real_end; i++) 
  {
    sum2 +=  v[i]*v[i]; 
  }

  return sum2; 

}

void TGraphAligned::getMeanAndRMS(Double_t * mean, Double_t * rms, Int_t istart, Int_t iend) const
{
  aligned_double_v v = GetY(); 
  __builtin_prefetch(v); 
  int N = GetN(); 
  int real_start = istart >= 0 ? TMath::Min(istart, N-1) : TMath::Max(0, N + istart)  ; 
  int real_end = iend >= 0 ? TMath::Min(iend, N-1) : TMath::Max(0, N + iend)  ; 

  double sum = 0; 
  double sum2 = 0; 


  //why the fuck is this not autovectorizing?!?? !?!?!? 
  //is it because hadd is too slow? 
  for (int i = real_start; i <= real_end; i++) 
  {
    sum +=  v[i]; 
    sum2 +=  v[i]*v[i]; 
  }

  double mn = sum/N; 
  if (mean) *mean = mn; 

  if (rms) 
  {
    double mn2 = mn*mn; 
    *rms = sqrt(sum2/N - mn2); 
  }
}


//TODO, this can be optimized mostly likely
double * TGraphAligned::getMoments(int N, double origin, double * moment) const
{
  if (!moment) moment = new double[N]; 


  for (int i = 0; i < N; i++ )
  {
    int n = i+1; 
    double sum = 0; 
    for (int j = 0; j < GetN(); j++) 
    {
      sum += TMath::Power(fX[j] - origin, n) * fY[j]; 
    }

    moment[i] = sum; 
  }


  return moment; 
}



void TGraphAligned::dBize(double mindB) 
{
  int n = GetN(); 

  for (int i = 0; i < n; i++)
  {
    fY[i] = TMath::Max(mindB, 10 * TMath::Log10(fY[i])); 
  }


}

Double_t TGraphAligned::pk2pk(Int_t nth_max, Int_t nth_min, Int_t * location_max, Int_t * location_min, Int_t istart, Int_t iend) const
{

  int start = istart < 0 ? GetN() + istart : istart; 
  int end = iend < 0 ? GetN() + iend : iend; 

  int max_index = -1; 
  int min_index = -1; 
  double max = 0; 
  double min = 0 ; 
  const double * y = GetY(); 

  
  //the easy case
  if (nth_max == 0 && nth_min == 0) 
  {
    for (int i = start; i <= end; i++) 
    {
      double val =y[i]; 
      if (max_index < 0 || val > max)
      {
        max = val; 
        max_index = i; 
      }

      if (min_index < 0 || val < min)
      {
        min = val; 
        min_index = i; 
      }
    }

    if (location_max) *location_max = max_index; 
    if (location_min) *location_min = min_index; 

  }

  else //probably not the optimal algorithm here. 
  {
    int N = end-start+1;
    double ycopy[N]; 
    memcpy(ycopy, y+start, N *sizeof(double)); 

    std::nth_element(ycopy, ycopy + N-1-nth_max, ycopy + N); 
    std::nth_element(ycopy, ycopy + nth_min, ycopy + N); 

    max = ycopy[N-1-nth_max]; 
    min = ycopy[nth_min]; 


    if (location_max || location_min)
    {
      bool need_max = location_max; 
      bool need_min = location_min;
      for (int i =start; i <=end; i++) 
      {
         if (location_max && y[i] == max) 
         {
           *location_max = i; 
           need_max = false;
         }

         if (location_min && y[i] == min) 
         {
           *location_min = i; 
           need_min = false;
         }

         if (!need_min && !need_max) break; 
      }
    }
  }

  return max-min; 
}

Double_t TGraphAligned::peakVal(Int_t * location, Int_t istart, Int_t iend, bool abs) const
{

  int start = istart < 0 ? GetN() + istart : istart; 
  int end = iend < 0 ? GetN() + iend : iend; 

  int max_index = -1; 
  double max = 0; 
  const double * y = GetY(); 

  for (int i = start; i <= end; i++) 
  {
    double val = abs ? fabs(y[i]) : y[i]; 
    if (max_index < 0 || val > max)
    {
      max = val; 
      max_index = i; 
    }
  }

  if (location) *location = max_index; 
  return max; 
}




