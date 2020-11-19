#include "TTree.h"
#include "TFile.h"
#include "FilterStrategy.h" 
#include "FilterOperation.h" 
#include "FilteredAnitaEvent.h" 


FilterStrategy::FilterStrategy(TFile * outfile) 
{
  started = false;
  f = 0; 
  if(outfile){
    attachFile(outfile);
  }
}


FilterStrategy::~FilterStrategy() 
{

  done(); 

  for (unsigned i = 0; i < operations.size(); i++) 
  {
    if (owns[i])
    {
      delete operations[i]; 
    }
  }

}

void FilterStrategy::attachFile(TFile* outfile){

  if(!outfile){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", outfile = " << outfile << std::endl;
    return;
  }
  
  if(f)
  {
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", FilterStrategy already has an outfile, ignoring passed outfile " << outfile << std::endl;
    return;
  }

  if(started)
  {
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to attach file after first call to FilterStrategy::process(FilteredAnitaEvent)" << std::endl;
    return;
  }

  f = outfile;

  // now we can make the trees with output branches since we have a file
  for(unsigned i=0; i < operations.size(); i++){
    FilterOperation* fo = operations[i];
    bool enable_output = enable_outputs[i];
	
    if (fo->nOutputs() && enable_output)     
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
  return;
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


void FilterStrategy::addOperation(FilterOperation* fo, bool enable_output, bool own) 
{
  if (started) 
  {
    //TODO: This really only needs to happen in the less common case that an operation has an output... but whatever. 
    fprintf(stderr,"Sorry, can't add operation to strategy after processing started!\n"); 
    return; 
  }

  operations.push_back(fo);
  enable_outputs.push_back(enable_output); 
  owns.push_back(own);

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


