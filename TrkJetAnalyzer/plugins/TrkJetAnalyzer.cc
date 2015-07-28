// -*- C++ -*-
//
// Package:    TrkJetAnalyzer
// Class:      TrkJetAnalyzer
// 
/**\class TrkJetAnalyzer TrkJetAnalyzer.cc TrkJetPhase2/TrkJetAnalyzer/plugins/TrkJetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Livia Soffi
//         Created:  Thu, 23 Jul 2015 12:55:55 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//---- for GenParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/Candidate/interface/Candidate.h"

//---- for GenJets
#include "DataFormats/JetReco/interface/GenJet.h" 

#include "TTree.h"
#include "TH1.h"
//
// class declaration
//


class TrkJetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrkJetAnalyzer(const edm::ParameterSet&);
      ~TrkJetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::InputTag GenJetCollection_;
  edm::InputTag GenParticlesCollection_;
  edm::InputTag TrackJetCollection_;
  edm::InputTag TrackCollection_;
  edm::InputTag VertexCollection_;
  edm::InputTag PUSummaryInfoCollection_;

  TTree* myTree_;

  //---Vertex and Pileup
 
  Int_t npu_;
  Float_t pu_zpos_[300];
  Float_t pu_sumpt_lowpt_[300];
  Float_t pu_sumpt_highpt_[300];
  Float_t pu_ntrks_lowpt_[300];
  Float_t pu_ntrks_highpt_[300];

  Float_t vxMC_;
  Float_t vyMC_;
  Float_t vzMC_;
  Int_t nvertex_;
  Float_t vx_[300];
  Float_t vy_[300];
  Float_t vz_[300];
  Float_t vntracks_[300];
  Float_t vchi2_[300];
  Float_t vndof_[300];




  //---- gen jets
  int ngenJet_;
  float genJetPt_[25];
  float genJetEta_[25];
  float genJetPhi_[25];
  float genJetMass_[25]; 
  float genJetEmE_[25]   ;
  float genJetHadrE_[25] ;
  float genJetInvE_[25]  ;
  float genJetAuxE_[25]  ;
  int genJetNconst_[25];

  //------- trk jets
  int ntrkJet_;
  float trkJetPt_[25];
  float trkJetEta_[25];
  float trkJetPhi_[25];
  float trkJetMass_[25]; 
 
  //------- tracks
  int ntrk_;
  float trkPt_[6000];
  float trkEta_[6000];
  float trkPhi_[6000];

  //----- gen particles
  int ngenCand_;
  float genCandPt_[30]   ;
  float genCandEta_[30]  ;
  float genCandPhi_[30]  ;
  float genCandMass_[30] ;
  float genCandStatus_[30];
  int genCandPdgId_[30];
  int genCandMothPdgId_[30];
  
   
  unsigned int minTracks_;
  TH1D *ntrkhisto;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrkJetAnalyzer::TrkJetAnalyzer(const edm::ParameterSet& iConfig):
  minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
{
   //now do what ever initialization is needed
  // I want to make a histogram of number of tracks in a<dwn event     
  edm::Service<TFileService> fs;

  GenJetCollection_       = iConfig.getParameter<edm::InputTag>("GenJetCollection");
  GenParticlesCollection_ = iConfig.getParameter<edm::InputTag>("GenParticlesCollection");
  TrackJetCollection_       = iConfig.getParameter<edm::InputTag>("TrackJetCollection");
  TrackCollection_       = iConfig.getParameter<edm::InputTag>("TrackCollection");
  VertexCollection_       = iConfig.getParameter<edm::InputTag>("VertexCollection");
  PUSummaryInfoCollection_ = iConfig.getParameter<edm::InputTag>("PUSummaryInfoCollection");

  myTree_ = fs -> make <TTree>("myTree","myTree");
  myTree_ -> Branch("ngenJet", &ngenJet_, "ngenJet/I");
  myTree_ -> Branch("genJetPt", genJetPt_, "genJetPt[30]/F");
  myTree_ -> Branch("genJetMass", genJetMass_, "genJet[30]/F");
  myTree_ -> Branch("genJetEta", genJetEta_, "genJetEta[30]/F");
  myTree_ -> Branch("genJetPhi", genJetPhi_, "genJetPhi[30]/F");
  myTree_ -> Branch("genJetEmE", genJetEmE_, "genJeEmEt[30]/F");
  myTree_ -> Branch("genJetHadrE", genJetHadrE_, "genJetHadrE[30]/F");
  myTree_ -> Branch("genJetInvE", genJetInvE_, "genJetInvE[30]/F");
  myTree_ -> Branch("genJetAuxE", genJetAuxE_, "genJetAuxE[30]/F");
  myTree_ -> Branch("genJetNconst", genJetNconst_, "genJetNconst[30]/I");
  
  myTree_ -> Branch("ntrkJet", &ntrkJet_, "ntrkJet/I");
  myTree_ -> Branch("trkJetPt", trkJetPt_, "trkJetPt[30]/F");
  myTree_ -> Branch("trkJetMass", trkJetMass_, "trkJet[30]/F");
  myTree_ -> Branch("trkJetEta", trkJetEta_, "trkJetEta[30]/F");
  myTree_ -> Branch("trkJetPhi", trkJetPhi_, "trkJetPhi[30]/F");
 
  myTree_ -> Branch("ntrk", &ntrk_, "ntrk/I");
  myTree_ -> Branch("trkPt", trkPt_, "trkPt[6000]/F");
  myTree_ -> Branch("trkEta", trkEta_, "trkEta[6000]/F");
  myTree_ -> Branch("trkPhi", trkPhi_, "trkPhi[6000]/F");
 

  myTree_ -> Branch("ngenCand",   &ngenCand_,     "ngenCand/I")   ;
  myTree_ -> Branch("genCandPt",   genCandPt_,    "genCandPt[30]/F")   ;
  myTree_ -> Branch("genCandEta" , genCandEta_,    "genCandEta[30]/F")  ;
  myTree_ -> Branch("genCandPhi",  genCandPhi_,    "genCandPhi[30]/F")  ;
  myTree_ -> Branch("genCandMass", genCandMass_,   "genCandMass[30]/F") ;
  myTree_ -> Branch("genCandStatus", genCandStatus_,  "genCandStatus[30]/F");
  myTree_ -> Branch("genCandPdgId", genCandPdgId_,  "genCandPdgId[30]/I");
  myTree_ -> Branch("genCandMothPdgId", genCandMothPdgId_,  "genCandMothPdgId[30]/I");
    
  myTree_->Branch("npu", &npu_, "npu/I");
  myTree_->Branch("pu_zpos", &pu_zpos_, "pu_zpos[npu]/F");
  myTree_->Branch("pu_sumpt_lowpt", &pu_sumpt_lowpt_, "pu_sumpt_lowpt[npu]/F");
  myTree_->Branch("pu_sumpt_highpt", &pu_sumpt_highpt_, "pu_sumpt_highpt[npu]/F");
  myTree_->Branch("pu_ntrks_lowpt", &pu_ntrks_lowpt_, "pu_ntrks_lowpt[npu]/F");
  myTree_->Branch("pu_ntrks_highpt", &pu_ntrks_highpt_, "pu_ntrks_highpt[npu]/F");


  myTree_->Branch("vxMC",&vxMC_,"vxMC/F");
  myTree_->Branch("vyMC",&vyMC_,"vyMC/F");
  myTree_->Branch("vzMC",&vzMC_,"vzMC/F");
  myTree_->Branch("nvertex",&nvertex_,"nvertex/I");
  myTree_->Branch("vx",&vx_,"vx[nvertex]/F");
  myTree_->Branch("vy",&vy_,"vy[nvertex]/F");
  myTree_->Branch("vz",&vz_,"vz[nvertex]/F");
  myTree_->Branch("vntracks",&vntracks_,"vntracks[nvertex]/F");
  myTree_->Branch("vchi2",&vchi2_,"vchi2[nvertex]/F");
  myTree_->Branch("vndof",&vndof_,"vndof[nvertex]/F");

  ntrkhisto = fs->make<TH1D>("tracks" , "Tracks" , 320 , 2500 , 10000 );

}


TrkJetAnalyzer::~TrkJetAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TrkJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   
   
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(GenParticlesCollection_,genParticles);
   
   edm::Handle<reco::GenJetCollection> genJet;
   iEvent.getByLabel(GenJetCollection_,genJet);
   
   edm::Handle<reco::TrackJetCollection> trkJet;
   iEvent.getByLabel(TrackJetCollection_,trkJet);
  
   edm::Handle<reco::TrackCollection> trk;
   iEvent.getByLabel(TrackCollection_,trk);
  
   edm::Handle<reco::VertexCollection> vtx;
   iEvent.getByLabel(VertexCollection_, vtx);

   edm::Handle<std::vector<PileupSummaryInfo> > pileup;
   iEvent.getByLabel(PUSummaryInfoCollection_, pileup);
      
  //save gen particles info
   ngenCand_ = 0;
   for(int i = 0;i<30;i++){
     genCandPt_[i] =-999.;
     genCandEta_[i]=-999.;
     genCandPhi_[i]=-999.;
     genCandMass_[i]=-999.;
     genCandStatus_[i]=-999.;
     genCandPdgId_[i]=-999;
     genCandMothPdgId_[i]=-999;
   }
   int k = 0;
   for (reco::GenParticleCollection::const_iterator genCandIter = genParticles->begin(); genCandIter != genParticles->end(); genCandIter++){
     if(k>30)continue;
     genCandPt_[k]= genCandIter->pt();
     genCandEta_[k]= genCandIter->eta();
     genCandPhi_[k]= genCandIter->phi();
     genCandMass_[k]= genCandIter->mass();
     genCandStatus_[k]= genCandIter->status();
     genCandPdgId_[k]= genCandIter->pdgId();
     // genCandMothPdgId_[k]= genCandIter->mother();//->pdgId();
     ngenCand_++;
     k++;

   }
   
   //save pu infos
   npu_ = 0;
 
   if( pileup.isValid() ) 
     {
       std::vector<PileupSummaryInfo>::const_iterator PVI;       
       for(PVI = pileup->begin(); PVI != pileup->end(); ++PVI) 
	 {
	   if(PVI->getBunchCrossing() != 0) 
	     continue;
	   npu_ = PVI->getPU_NumInteractions();
	   int sv = PVI->getPU_zpositions().size() < 50 ? PVI->getPU_zpositions().size() : 300;
	   for (int iPU=0;iPU<sv;++iPU)
	     {
	       pu_zpos_[iPU]=PVI->getPU_zpositions()[iPU];
	       pu_sumpt_lowpt_[iPU]=PVI->getPU_sumpT_lowpT()[iPU];
	       pu_sumpt_highpt_[iPU]=PVI->getPU_sumpT_highpT()[iPU];
	       pu_ntrks_lowpt_[iPU]=PVI->getPU_ntrks_lowpT()[iPU];
	       pu_ntrks_highpt_[iPU]=PVI->getPU_ntrks_highpT()[iPU];
	     }
	 }
     }
   
   
  // Get the primary vertex coordinates
   nvertex_=0;
   for (reco::VertexCollection::const_iterator it = vtx->begin(); 
	it != vtx->end(); ++it) {
     
     vx_[nvertex_] = (it->isValid()) ? it->x() : 999.;
     vy_[nvertex_] = (it->isValid()) ? it->y() : 999.;
     vz_[nvertex_] = (it->isValid()) ? it->z() : 999.;
     
     vntracks_[nvertex_] = (it->isValid()) ? it->tracksSize() : 0;
     vchi2_[nvertex_] = (it->isValid()) ? it->normalizedChi2() : 100.;
     vndof_[nvertex_] = (it->isValid()) ? it->ndof() : 0.;
     
     nvertex_++;
   }
   
   //save gen jet informations
   ngenJet_ = 0;
   for(int i = 0;i<25;i++){
     genJetPt_[i] =-999.;
     genJetEta_[i]=-999.;
     genJetPhi_[i]=-999.;
     genJetMass_[i]=-999.;
     genJetEmE_[i]=-999.;
     genJetHadrE_[i]=-999.;
     genJetInvE_[i]=-999.;
     genJetAuxE_[i]=-999.;
     genJetNconst_[i]=-999;
 }
   int j=0;
   for (reco::GenJetCollection::const_iterator genJetIter=genJet->begin(); genJetIter!=genJet->end(); genJetIter++){
     if(genJetIter->pt()<5)continue; 
     genJetPt_[j]=genJetIter->pt();
     genJetEta_[j]=genJetIter->eta();
     genJetPhi_[j]=genJetIter->phi();
     genJetMass_[j]=genJetIter->mass();
     genJetEmE_[j]=genJetIter->emEnergy();
     genJetHadrE_[j]  =genJetIter->hadEnergy();
     genJetInvE_[j]=genJetIter->invisibleEnergy();
     genJetAuxE_[j]=genJetIter->auxiliaryEnergy();
     genJetNconst_[j]=genJetIter->nConstituents();
     ngenJet_++;
     j++;
   }
   
   //save trk jet infos
   ntrkJet_=0;
   for(int l = 0;l<25;l++){
     trkJetPt_[l] =-999.;
     trkJetEta_[l]=-999.;
     trkJetPhi_[l]=-999.;
     trkJetMass_[l]=-999.;
 }
   int z=0;
   for (reco::TrackJetCollection::const_iterator trkJetIter=trkJet->begin(); trkJetIter!=trkJet->end(); trkJetIter++){
     if(trkJetIter->pt()<5)continue; 
     trkJetPt_[z]=trkJetIter->pt();
     trkJetEta_[z]=trkJetIter->eta();
     trkJetPhi_[z]=trkJetIter->phi();
     trkJetMass_[z]=trkJetIter->mass();
     z++;
     ntrkJet_++;
   }
  
   //save trk  infos
   ntrk_=0;
   for(int ll = 0;ll<6000;ll++){
     trkPt_[ll] =-999.;
     trkEta_[ll]=-999.;
     trkPhi_[ll]=-999.;
   }
   int zz=0;
   for (reco::TrackCollection::const_iterator trkIter=trk->begin(); trkIter!=trk->end(); trkIter++){
     //if(trkIter->pt()<5)continue; 
     trkPt_[zz]=trkIter->pt();
     trkEta_[zz]=trkIter->eta();
     trkPhi_[zz]=trkIter->phi();
     zz++;
     ntrk_++;
   }
  
  
   
   myTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrkJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrkJetAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TrkJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TrkJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TrkJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TrkJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrkJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkJetAnalyzer);
