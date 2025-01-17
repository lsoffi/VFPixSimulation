
// last commit by $Id: analysis.cc,v 1.1 2013/01/31 15:32:00 soffi Exp $
//
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TAxis.h>

#include "Analysis.h"

using namespace std;



int main(int argc, char* argv[]) {


      //================ Creating chain 
  

  TFile* fin= new TFile("../TrkJetAnalyzer/histoOUTPUT.root");
  TChain* chain =new TChain("myTree");

  chain->Add("../TrkJetAnalyzer/histoOUTPUT.root/TrkJetPhase2/myTree");


  // std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run analysis
  Analysis tree( chain );
  tree.Loop();
  
  delete fin;
  delete chain;
  return 0;
}
