#ifndef _ANITA_TMVA_H
#define _ANITA_TMVA_H


#include <vector>
#include <TMVA/Reader.h>
#include <TMVA/Factory.h>


class TTree; 
class TFile; 

/** 
 *
 * TMVA helpers for ANITA analysis 
 *
 *  The basic problem is that AnitaEventSummary has arrays of things, but TMVA
 *  doesn't support that properly (even though... it could in a way similar to
 *  TTree::Scan) 
 *
 *  See macros/TMVAExample.C for an example analysis. 
 *
 *
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 *
 */
 
namespace AnitaTMVA
{

  /** A variable to be used for multivariable analysis */
  struct MVAVar
  {
    MVAVar(const char * expr, const char * name, char type='F', bool spectator=false)  :  expression(expr), name(name), type(type), spectator(spectator)  {;} 
    const char * expression;  /* A valid formula for this variable (made with TTree::Draw) */
    const char * name;  /* a name you assign to this variable  (no spaces, will be used for branch name( */ 
    char type; /*F' or 'I'. TMVA converts everything to float internally, so doubles will be converted to floats. Sorry.  */
    bool spectator; 
  }; 


  /** A set of variables for multivariable analysis */ 
  struct MVAVarSet
  {

    /* empty var set */ 
    MVAVarSet() { ;} 

    /* make from arrays */ 
    MVAVarSet( int nvars,  const char * varNames[], const char * varExpr[], const char * type = 0, const bool * spectator = 0); //defaults to 'F'  

    /* Use this TMVA::DataLoader (or a TMVA::Factory if using an old TMVA) */ 
    template <typename T> 
    void setUpData(T *t) const
    { 
      for (unsigned i = 0; i < vars.size(); i++) 
      {
        if (!vars[i].spectator)
          t->AddVariable(vars[i].name, vars[i].type); 
        else
          t->AddSpectator(vars[i].name, vars[i].type); 
      }
    }

    void add(const MVAVar& var) { vars.push_back(var); }

    int N() const { return (int) vars.size(); } 

    const MVAVar & at(int i) const { return vars[i]; } 

    std::vector<MVAVar> vars; 
  }; 


  /** Creates a tree that may be used by TMVA. Stores in in file outfile, which must exist. 
   *
   * This tree can then be used for training or reading using the varName you provide.  Note that internally this uses TTree::Draw, so if your cuts suck, you might
   * end up with an insane amount of memory usage, but since TMVA loads everything into memory anyway, I don't think you lose much. This is necessary because TMVA's 
   * methods with selection don't properly handle cuts on array variables. Maybe in the future they'll fix this and this won't be necessary. 
   *
   * The tree will also have the iteration and entry number in it. 
   *
   * @param intree The input tree  
   * @param outfile The file to write the output tree to
   * @param out_tree_hname the name of the new tree
   * @param varsAn MVAVarSet describing the variables to use 
   * @param selection The selection cuts. 
   * @return pointer to the created tree
   *
   */


  TTree* makeTMVATree(TTree * intree, 
                   TFile * outfile,
                   const char * out_tree_name, 
                   const MVAVarSet & vars,
                    const char * selection
                   ); 


  /** Adds the TMVA classifier result to a tree.  This tree should have been
   * generated with makeTMVATree. A branch will be added to the tree with the
   * same name as the method 
   *
   * This takes care of generating the reader and such. 
   *
   * @param tree the tree to evaluate. Should have been generated with maketMVATree
   * @param vars The set of variables that are being evaluated. All branches must exist in tree and 
   *             the TMVA method must have been generated with it. 
   * @param branch_name The name of the evaluation branch to be added to the tree
   * @param tmva_weights_file The weights.xml file that you want to read in. 
   * @param aux aux parameter to EvaluateMVA
   * @reuturn 0 on success;
   * */

  int evaluateTMVA(TTree * tree, const MVAVarSet & vars,
                    const char * branch_name,
                    const char * tmva_weights_xml_file,
                    double aux = 0) ; 


    

}
#endif
