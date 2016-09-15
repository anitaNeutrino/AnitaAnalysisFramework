#include "TTree.h"
#include "TFile.h"
#include "FilterStrategy.h" 
#include "FilterOperation.h" 
#include "FilteredAnitaEvent.h" 


FilterStrategy::FilterStrategy(TFile * outfile) 
{
  f = outfile; 
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
    if (f && operations[op]->nOutputs() > 0 && enable_outputs[op])
    {
      f->cd(); 
      for (unsigned j = 0; j < operations[op]->nOutputs(); j++)
      {
        operations[op]->fillOutput(j,&outputStore[tree_index][j][0]); 
      }
      trees[tree_index]->Fill(); 
      tree_index++; 
    }
  }
}


void FilterStrategy::addOperation(FilterOperation* fo, bool enable_output) 
{
  if (started) 
  {
    fprintf(stderr,"Sorry, can't add operation to strategy after processing started!\n"); 
    return; 
  }

  operations.push_back(fo);
  enable_outputs.push_back(enable_output); 

  if (f && fo->nOutputs() && enable_output) 
  {
    const char * tag = fo->tag(); 
    std::string stag(tag); 
    int i = used_ids.count(stag); 
    used_ids.insert(stag); 
    f->cd(); 


    TString tagstr = i ? TString::Format("%s%d",tag,i) : TString(tag); 
    TTree * tree = new TTree(tagstr.Data(), fo->description()); 
    trees.push_back(tree); 
    std::vector<std::vector<double> > store(fo->nOutputs()); 
    size_t which_store = outputStore.size();
    outputStore.push_back(store); 

    for (unsigned j = 0; j < fo->nOutputs(); j++) 
    {
      outputStore[which_store][j].insert(outputStore[which_store][j].end(),fo->outputLength(j),0); 
      tree->Branch(fo->outputName(j), &outputStore[which_store][j][0], TString::Format("%s[%d]/D",fo->outputName(j), fo->outputLength(j))); 
    }
  }
}


