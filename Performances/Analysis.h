//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 28 13:40:51 2015 by ROOT version 5.34/07
// from TTree myTree/myTree
// found on file: ../TrkJetAnalyzer/histoOUTPUT.root
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ngenJet;
   Float_t         genJetPt[30];
   Float_t         genJetMass[30];
   Float_t         genJetEta[30];
   Float_t         genJetPhi[30];
   Float_t         genJetEmE[30];
   Float_t         genJetHadrE[30];
   Float_t         genJetInvE[30];
   Float_t         genJetAuxE[30];
   Int_t           genJetNconst[30];
   Int_t           ntrkJet;
   Float_t         trkJetPt[30];
   Float_t         trkJetMass[30];
   Float_t         trkJetEta[30];
   Float_t         trkJetPhi[30];
   Int_t           ntrk;
   Float_t         trkPt[6000];
   Float_t         trkEta[6000];
   Float_t         trkPhi[6000];
   Int_t           ngenCand;
   Float_t         genCandPt[30];
   Float_t         genCandEta[30];
   Float_t         genCandPhi[30];
   Float_t         genCandMass[30];
   Float_t         genCandStatus[30];
   Int_t           genCandPdgId[30];
   Int_t           genCandMothPdgId[30];
   Int_t           npu;
   Float_t         pu_zpos[182];   //[npu]
   Float_t         pu_sumpt_lowpt[182];   //[npu]
   Float_t         pu_sumpt_highpt[182];   //[npu]
   Float_t         pu_ntrks_lowpt[182];   //[npu]
   Float_t         pu_ntrks_highpt[182];   //[npu]
   Float_t         vxMC;
   Float_t         vyMC;
   Float_t         vzMC;
   Int_t           nvertex;
   Float_t         vx[113];   //[nvertex]
   Float_t         vy[113];   //[nvertex]
   Float_t         vz[113];   //[nvertex]
   Float_t         vntracks[113];   //[nvertex]
   Float_t         vchi2[113];   //[nvertex]
   Float_t         vndof[113];   //[nvertex]

   // List of branches
   TBranch        *b_ngenJet;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJet;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJeEmEt;   //!
   TBranch        *b_genJetHadrE;   //!
   TBranch        *b_genJetInvE;   //!
   TBranch        *b_genJetAuxE;   //!
   TBranch        *b_genJetNconst;   //!
   TBranch        *b_ntrkJet;   //!
   TBranch        *b_trkJetPt;   //!
   TBranch        *b_trkJet;   //!
   TBranch        *b_trkJetEta;   //!
   TBranch        *b_trkJetPhi;   //!
   TBranch        *b_ntrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_ngenCand;   //!
   TBranch        *b_genCandPt;   //!
   TBranch        *b_genCandEta;   //!
   TBranch        *b_genCandPhi;   //!
   TBranch        *b_genCandMass;   //!
   TBranch        *b_genCandStatus;   //!
   TBranch        *b_genCandPdgId;   //!
   TBranch        *b_genCandMothPdgId;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_pu_zpos;   //!
   TBranch        *b_pu_sumpt_lowpt;   //!
   TBranch        *b_pu_sumpt_highpt;   //!
   TBranch        *b_pu_ntrks_lowpt;   //!
   TBranch        *b_pu_ntrks_highpt;   //!
   TBranch        *b_vxMC;   //!
   TBranch        *b_vyMC;   //!
   TBranch        *b_vzMC;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vntracks;   //!
   TBranch        *b_vchi2;   //!
   TBranch        *b_vndof;   //!

   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../TrkJetAnalyzer/histoOUTPUT.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../TrkJetAnalyzer/histoOUTPUT.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../TrkJetAnalyzer/histoOUTPUT.root:/TrkJetPhase2");
      dir->GetObject("myTree",tree);

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ngenJet", &ngenJet, &b_ngenJet);
   fChain->SetBranchAddress("genJetPt", genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetMass", genJetMass, &b_genJet);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetEmE", genJetEmE, &b_genJeEmEt);
   fChain->SetBranchAddress("genJetHadrE", genJetHadrE, &b_genJetHadrE);
   fChain->SetBranchAddress("genJetInvE", genJetInvE, &b_genJetInvE);
   fChain->SetBranchAddress("genJetAuxE", genJetAuxE, &b_genJetAuxE);
   fChain->SetBranchAddress("genJetNconst", genJetNconst, &b_genJetNconst);
   fChain->SetBranchAddress("ntrkJet", &ntrkJet, &b_ntrkJet);
   fChain->SetBranchAddress("trkJetPt", trkJetPt, &b_trkJetPt);
   fChain->SetBranchAddress("trkJetMass", trkJetMass, &b_trkJet);
   fChain->SetBranchAddress("trkJetEta", trkJetEta, &b_trkJetEta);
   fChain->SetBranchAddress("trkJetPhi", trkJetPhi, &b_trkJetPhi);
   fChain->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("ngenCand", &ngenCand, &b_ngenCand);
   fChain->SetBranchAddress("genCandPt", genCandPt, &b_genCandPt);
   fChain->SetBranchAddress("genCandEta", genCandEta, &b_genCandEta);
   fChain->SetBranchAddress("genCandPhi", genCandPhi, &b_genCandPhi);
   fChain->SetBranchAddress("genCandMass", genCandMass, &b_genCandMass);
   fChain->SetBranchAddress("genCandStatus", genCandStatus, &b_genCandStatus);
   fChain->SetBranchAddress("genCandPdgId", genCandPdgId, &b_genCandPdgId);
   fChain->SetBranchAddress("genCandMothPdgId", genCandMothPdgId, &b_genCandMothPdgId);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("pu_zpos", pu_zpos, &b_pu_zpos);
   fChain->SetBranchAddress("pu_sumpt_lowpt", pu_sumpt_lowpt, &b_pu_sumpt_lowpt);
   fChain->SetBranchAddress("pu_sumpt_highpt", pu_sumpt_highpt, &b_pu_sumpt_highpt);
   fChain->SetBranchAddress("pu_ntrks_lowpt", pu_ntrks_lowpt, &b_pu_ntrks_lowpt);
   fChain->SetBranchAddress("pu_ntrks_highpt", pu_ntrks_highpt, &b_pu_ntrks_highpt);
   fChain->SetBranchAddress("vxMC", &vxMC, &b_vxMC);
   fChain->SetBranchAddress("vyMC", &vyMC, &b_vyMC);
   fChain->SetBranchAddress("vzMC", &vzMC, &b_vzMC);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vntracks", vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", vchi2, &b_vchi2);
   fChain->SetBranchAddress("vndof", vndof, &b_vndof);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
