#include "TGraphAligned.h" 

#include <stdlib.h>
#include <malloc.h>
#include "TH1.h"
#include "TList.h"


ClassImp(TGraphAligned); 


TGraphAligned::TGraphAligned()
  :TGraph()
{
  //don't think we need to do anything in this case. 
}

TGraphAligned::TGraphAligned(Int_t n)
{
  SetTitle("Graph"); 
  SetName("Graph"); 
  fNpoints =n; 
  if (!CtorAllocate()) return; 

  FillZero(0, fNpoints); 
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
      } else {
        fX = (Double_t*) memalign(ALIGNMENT, fNpoints * sizeof(Double_t)); 
        fY = (Double_t*) memalign(ALIGNMENT, fNpoints * sizeof(Double_t)); 
      }

      Int_t n = gr.GetN() * sizeof(Double_t);
      if (n > 0) {
         memcpy(fX, gr.fX, n);
         memcpy(fY, gr.fY, n);
      }
   }
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
        fX = (Double_t*) memalign(ALIGNMENT, fNpoints * sizeof(Double_t)); 
        fY = (Double_t*) memalign(ALIGNMENT, fNpoints * sizeof(Double_t)); 
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
      Double_t * new_fX = (Double_t*) memalign(ALIGNMENT, fMaxSize * sizeof(Double_t)); 
      memcpy(new_fX, fX, gr.GetN() * sizeof(Double_t));
      delete [] fX; 
      fX = new_fX; 

      Double_t * new_fY = (Double_t*) memalign(ALIGNMENT, fMaxSize * sizeof(Double_t)); 
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
       newarrays[i] = (Double_t*) memalign(ALIGNMENT, arraySize * sizeof(Double_t)); 
    }
  }

  fMaxSize = arraySize;
  return newarrays;
}
