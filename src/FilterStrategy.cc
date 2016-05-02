#include "TTree.h"
#include "TFile.h"
#include "FilterStrategy.h" 
#include "FilterOperation.h" 
#include "FilteredAnitaEvent.h" 


FilterStrategy::FilterStrategy(TFile * outfile) 
{
  if (outfile) 
  {
    f =  outfile; 
  }
  started = false; 
}

void FilterStrategy::done() 
{
  if (f) 
  {
    f->cd(); 

    for (size_t i = 0; i < trees.size(); i++) 
    {
      trees[i]->Write(); 
    }

    f->Flush(); 
  }


  for (size_t i = 0; i < trees.size(); i++) 
  {
     delete trees[i]; 
  }
}

void FilterStrategy::process(FilteredAnitaEvent * ev) 
{
  started = true; 
  int tree_index = 0; 


  for (size_t op = 0; op < operations.size(); op++) 
  {
    //Do the thing! 
    operations[op]->process(ev); 


    // save the stage if are asked to  it's not the last... the last is saved anyway
    if (ev->keep_all_stages && op < operations.size() - 1)
    {
      ev->saveStage(operations.size()-1); 
    }

    // fill the output trees 
    if (f && operations[op]->nOutputs() > 0)
    {
      f->cd(); 
      operations[op]->fillOutputs(&outputStore[tree_index][0]); 
      trees[tree_index]->Fill(); 
      tree_index++; 
    }
  }
}


void FilterStrategy::addOperation(FilterOperation* fo) 
{
  if (started) 
  {
    fprintf(stderr,"Sorry, can't add operation to strategy after processing started!\n"); 
    return; 
  }

  operations.push_back(fo);

  if (f && fo->nOutputs()) 
  {
    const char * tag = fo->tag(); 
    std::string stag(tag); 
    int i = used_ids.count(stag); 
    used_ids.insert(stag); 
    f->cd(); 


    TString tagstr = i ? TString::Format("%s%d",tag,i) : TString(tag); 
    TTree * tree = new TTree(tagstr.Data(), fo->description()); 
    trees.push_back(tree); 
    std::vector<double> store(fo->nOutputs()); 
    outputStore.push_back(store); 

    for (unsigned j = 0; j < fo->nOutputs(); j++) 
    {
      tree->Branch(fo->outputName(j), &outputStore[outputStore.size()-1][j]); 
    }
  }
}


