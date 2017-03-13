#include "AnitaTMVA.h" 

/*  Macro to demonstrate usage of TMVA with AnitaEventSummary's
 *  
 *  Note that the TMVA interface changed completely between ROOT 6.06 and ROOT 6.08. This macro
 *  defaults to the ROOT 6.08 interface, but it's easy to modify for ROOT 6.06. 
 *
 *  Please don't use this for real analysis, it's just meant to chose how to do things. In this case, cut_signal selects for a WAIS pulser and cut_bg selects for not a wais pulser  
 *
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 */ 


/* These are the names of our branches */ 
 
const char * vars[] = { "mapPeak","hilbertPeak","peakTime" }; 

/* These are the expressions for our branches */ 
const char * expr[] = {"peak.value[][]", "coherent.peakHilbert[][]", "coherent.peakTime[][]"}; 

/* These is our selection cuts.
 *
 * This is just an example, so I want to show how to both train TMVA and use it to apply cuts, so I'll use odd events for one and even for the other 
 *
 */
const char * cut_signal = "(Entry$ %2) == 1  && peak.value[][] > 0 && flags.isRF==1 && peak[][].theta < 60 && peak[][].theta > -60 && flags.pulser == 1";
const char * cut_bg = "(Entry$ % 2) ==1  && peak.value[][] > 0 && flags.isRF==1 && peak[][].theta < 60 && peak[][].theta > -60 && flags.pulser == 0";
const char * cut_eval = "(Entry$ %2 ) == 0 && peak.value[][] > 0 && flags.isRF==1 && peak[][].theta < 60 && peak[][].theta > -60 ";


/* filename should be a ROOT file containg tree with AnitaEventSummary treename */ 
void TMVAExample(const char * filename, const char * treename)
{

  /* This is just a lazy way of loading the tree in this case, but you could also add multiple files trivially */ 
  TChain c(treename); 
  c.Add(filename); 

  /* Create our variable set */ 
  AnitaTMVA::MVAVarSet varset(3,vars,expr); 

  /* We'll do everything in the same file here. */ 
  TFile f("tmva.root","RECREATE"); 
  TTree *sigtree = AnitaTMVA::makeTMVATree(&c, &f, "signal_in", varset, cut_signal); 
  TTree *bgtree = AnitaTMVA::makeTMVATree(&c, &f, "bg_in", varset, cut_bg); 
  TTree *evaltree = AnitaTMVA::makeTMVATree(&c, &f, "eval_in", varset, cut_eval); 

  /** Create the TMVA factory */ 
  TMVA::Factory factory("tmva_example", &f); 


  /** If you're using ROOT < 6.08, comment out everything with the data_loader and uncomment the following lines out */ 
  // varset.setUpData(&factory); 
  // factory.AddSignalTree(sigtree); 
  // factory.AddBackgroundTree(bgtree); 
  // factory.BookMethod(TMVA::Types::kFisher,"Fisher"); // book a Fisher discriminant 


  TMVA::DataLoader loader; 
  varset.setUpData(&loader); 
  loader.AddSignalTree(sigtree); 
  loader.AddBackgroundTree(bgtree); 
  factory.BookMethod(&loader, TMVA::Types::kFisher,"Fisher"); // book a Fisher discriminant 


  //Tell TMVA to train, test and evaluate
  factory.TrainAllMethods(); 
  factory.TestAllMethods(); 
  factory.EvaluateAllMethods(); 


  // now let's go ahead and tag our eval tree with the output. 
  // not sure why this is the default path on my system, maybe it's different on yours? 
  evaluateTMVA( evaltree, varset, "Fisher",  "default/weights/tmva_example_Fisher.weights.xml"); 

  f.Write(); 
}

