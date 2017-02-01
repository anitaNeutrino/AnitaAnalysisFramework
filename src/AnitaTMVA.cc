#include "AnitaTMVA.h" 
#include "TTree.h" 
#include "TFile.h" 
#include <string>
#include <sstream>



union idiocy
{
  int i;
  float f;
}; 


AnitaTMVA::MVAVarSet::MVAVarSet(int n, const char * names[], const char * expr[], const char * type, const bool * spectator) 
{
  for (int i = 0; i < n; i++) 
  {
    add(MVAVar(expr[i], names[i], type ? type[i] :'F', spectator ? spectator[i] : false)); 
  }
}


TTree* AnitaTMVA::makeTMVATree(TTree * in, TFile * outf, const char * tree_name, const AnitaTMVA::MVAVarSet & vars , const char * cut) 
{

  std::stringstream drawstr; 

  outf->cd(); 
  TTree * out = new TTree(tree_name,tree_name); 


  std::vector<idiocy> mem(vars.N()+2); 

  for (int i = 0; i < vars.N(); i++) 
  {
    drawstr << vars.at(i).expression << ":"; 

    char type = vars.at(i).type; 
    switch(type)
    {
      case 'I': 
        out->Branch(vars.at(i).name, &mem[i].i); break; 
      case 'F': 
      default:
        out->Branch(vars.at(i).name, &mem[i].f);
    }
  }

  out->Branch("entry", &mem[vars.N()].i); 
  out->Branch("iteration", &mem[vars.N()+1].i); 
  drawstr << "Entry$:Iteration$"; 

//  printf("%s\n",drawstr.str().c_str()); 

  //now the real work happens... 
  int Nout = in->Draw(drawstr.str().c_str(), cut,"goff"); 
 

  //This is madness. I should just use a custom selector but I'm too lazy for that right now. 
  for (int j =0; j < Nout; j++) 
  {

    for (int i = 0; i < vars.N()+2; i++)
    {
      char type = i < vars.N() ? vars.at(i).type : 'I'; 
      switch(type) 
      {
        case 'I':
          mem[i].i = in->GetVal(i)[j]; break;
        case 'F': 
        default: 
          mem[i].f = in->GetVal(i)[j]; 
      }
    }

    out->Fill(); 
  }

  return out; 
}


int AnitaTMVA::evaluateTMVA(TTree * tree, const AnitaTMVA::MVAVarSet & vars, const char * branch_name, const char * weights_file, double aux) 
{
  std::vector<idiocy> mem(vars.N()); 
  TMVA::Reader reader; 


  for (int i = 0; i < vars.N(); i++) 
  {
    if(! tree->FindBranch(vars.at(i).name))
    {
      fprintf(stderr,"%s not found in tree... did you create this tree with the same variable set? Aborting evaluateTMVA.", vars.at(i).name); 
      return 1; 
    }

    if (vars.at(i).spectator) continue; 

    switch(vars.at(i).type)
    {
      case 'I':
         tree->SetBranchAddress(vars.at(i).name, &mem[i].i); 
         reader.AddVariable(vars.at(i).name, &mem[i].i); 
      case 'F':
      default:
         tree->SetBranchAddress(vars.at(i).name, &mem[i].f); 
         reader.AddVariable(vars.at(i).name, &mem[i].f); 
         break;
    }
  }

  float value = 0; 
  reader.BookMVA(branch_name,weights_file); 
  TBranch * b = tree->Branch( branch_name ,&value); 

  for (int j = 0; j < tree->GetEntries(); j++) 
  {
    tree->GetEntry(j); 
    value = reader.EvaluateMVA(branch_name,aux); 
    b->Fill(); 
  }



  return 0; 

}





