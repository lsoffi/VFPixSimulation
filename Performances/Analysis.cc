#define Analysis_cxx
#include "Analysis.h"
#include"TPaveText.h"
#include "TChain.h"
#include "TH1F.h"
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"                                                                                                                                                         
#include "TLegend.h"                                                                                                                                                        
#include "TPaveText.h"                                                                                                                                                      
#include "TColor.h"                                                                                                                                                         
#include "TLatex.h"                                                                                                                                                         
#include "TLorentzVector.h"

//#include "../CMS_lumi.C"

void Analysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   bool isMatched[300];
   Double_t trkJetDR[300];
   Double_t genJetAssocId[300];

   TH1F* h_ngenJet = new TH1F("h_ngenJet","h_ngenJet", 200,-0.5,200.5);
   TH1F* h_ntrkJet = new TH1F("h_ntrkJet","h_ntrkJet", 200,-0.5,200.5);
   TH1F* h_ntrk = new TH1F("h_ntrk","h_ntrk", 200,-0.5,200.5);

   TH1F* h_genJetPt = new TH1F("h_genJetPt","h_genJetPt", 160,0,1200);
   TH1F* h_genJetEta = new TH1F("h_genJetEta","h_genJetEta", 160,-10,10);
   TH1F* h_genJetPhi = new TH1F("h_genJetPhi","h_genJetPhi", 160,-6,6);
   TH1F* h_genJetMass = new TH1F("h_genJetMass","h_genJetMass", 160,0,1200);
  
   TH1F* h_trkJetPt = new TH1F("h_trkJetPt","h_trkJetPt", 160,0,1200);
   TH1F* h_trkJetEta = new TH1F("h_trkJetEta","h_trkJetEta", 160,-10,10);
   TH1F* h_trkJetPhi = new TH1F("h_trkJetPhi","h_trkJetPhi", 160,-6,6);
   TH1F* h_trkJetMass = new TH1F("h_trkJetMass","h_trkJetMass", 160,0,1200);
   TH1F* h_trkJetDR = new TH1F("h_trkJetDR","h_trkJetDR", 160,0,10);

   TH2F* h2_pt_eta = new TH2F("h2_pt_eta","h2_pt_eta",160,0,10,160, 0., 1200);
   TH2F* h2_phi_eta = new TH2F("h2_phi_eta","h2_phi_eta",160,0,10,160, -6., 6);
   TH2F* h2_resp_eta = new TH2F("h2_resp_eta","h2_resp_eta",160,0,10,60, 0., 10);
   

   TH1F* h_trkPt = new TH1F("h_trkPt","h_trkPt", 160,0,2000);
   TH1F* h_trkEta = new TH1F("h_trkEta","h_trkEta", 160,-6,6);
   TH1F* h_trkPhi = new TH1F("h_trkPhi","h_trkPhi", 160,-6,6);
   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%100==0)std::cout<<jentry<<std::endl;
      
      //match gen and trk jets
      for(int trk = 0; trk< ntrkJet;trk++){
	if(abs(trkJetEta[trk])>5)continue;
	std::cout<<"trk: "<<trk<<std::endl;
	double dRMin = 999.;
	int idGen=299;

	isMatched[trk]=false;
	trkJetDR[trk]=999.;

	for(int gen = 0; gen< ngenJet;gen++){
	  std::cout<<"gen: "<<gen<<std::endl;
	  double dRJet = 999.;

	  TLorentzVector* trkP4 = new TLorentzVector();
	  trkP4->SetPtEtaPhiM(trkJetPt[trk],trkJetEta[trk],trkJetPhi[trk],trkJetMass[trk]);
	  TLorentzVector* genP4 = new TLorentzVector();
	  genP4->SetPtEtaPhiM(genJetPt[gen],genJetEta[gen],genJetPhi[gen],genJetMass[gen]);
	  if( genP4->Pt() < 10. || fabs(genP4->Eta()) > 10. ) continue;
	  dRJet = trkP4->DeltaR(*genP4);
	  std::cout<<"dRJet: " <<dRJet<<"  id: "<<gen<<std::endl;
	  if (dRJet<dRMin/*&&fabs((trkP4->Pt()-genP4->Pt())/trkP4->Pt()) < 0.5*/ ){	    
	    dRMin=dRJet;
	    idGen = gen; 
	  }
	}
	trkJetDR[trk]=dRMin;
	genJetAssocId[trk]=idGen;
	std::cout<<"----------->chosen: "<<idGen<<"   dR: "<<trkJetDR[trk]<<std::endl;
	      	
      }
      
      h_ntrk->Fill(ntrk);
      std::cout<<"fill histos for trk"<<std::endl; 
      for(int i = 0;i<ntrk;i++){
	h_trkPt->Fill(trkPt[i]);
	h_trkEta->Fill(trkEta[i]);
	h_trkPhi->Fill(trkPhi[i]);
   }

      
      h_ntrkJet->Fill(ntrkJet);
      std::cout<<"fill histos for trk Jet"<<std::endl;
      for(int ii = 0;ii<ntrkJet;ii++){
	h_trkJetPt->Fill(trkJetPt[ii]);
	h_trkJetEta->Fill(trkJetEta[ii]);
	h_trkJetPhi->Fill(trkJetPhi[ii]);
	h_trkJetMass->Fill(trkJetMass[ii]);
	h_trkJetDR->Fill(trkJetDR[ii]);
	h2_pt_eta->Fill(fabs(trkJetEta[ii]),trkJetPt[ii]);
	h2_phi_eta->Fill(fabs(trkJetEta[ii]),trkJetPhi[ii]);
	int genId = genJetAssocId[ii];
	if(genId<30){	  
	  h2_resp_eta->Fill(fabs(trkJetEta[ii]),trkJetPt[ii]/genJetPt[genId]);
	}
      }
       
      h_ngenJet->Fill(ngenJet);
       std::cout<<"fill histos for gen Jet"<<std::endl;
      for(int i = 0;i<ngenJet;i++){
	h_genJetPt->Fill(genJetPt[i]);
	h_genJetEta->Fill(genJetEta[i]);
	h_genJetPhi->Fill(genJetPhi[i]);
	h_genJetMass->Fill(genJetMass[i]);
      }

    
   }
   gStyle->SetOptStat(0);
   TCanvas* c1 = new TCanvas("c1", "c1", 1);
   c1->cd();
   //CMS_lumi(c1,true,0);
 
   //trk pt
   h_trkPt->SetLineColor(kAzure+7);
   h_trkPt->Draw("hist"); 
   h_trkPt->GetYaxis()->SetTitle("Events");     
   h_trkPt->GetXaxis()->SetTitle("Tracks p_{T} [GeV]");     
   c1->SetLogy(0); 
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Trakcs_Pt.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Trakcs_Pt_LOG.png");
   
   //trk eta
   h_trkEta->SetLineColor(kAzure+7);
   h_trkEta->Draw("hist"); 
   h_trkEta->GetYaxis()->SetTitle("Events");     
   h_trkEta->GetXaxis()->SetTitle("Tracks #eta");     
   c1->SetLogy(0); 
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Trakcs_Eta.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Trakcs_Eta_LOG.png");
   
   //trk phi
   h_trkPhi->SetLineColor(kAzure+7);
   h_trkPhi->Draw("hist"); 
   h_trkPhi->GetYaxis()->SetTitle("Events");     
   h_trkPhi->GetXaxis()->SetTitle("Tracks #eta");     
   c1->SetLogy(0);   
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Trakcs_Phi.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Trakcs_Phi_LOG.png");
   


 
   //jet n
   h_ntrkJet->SetLineColor(kAzure+7);
   h_ngenJet->SetLineColor(kMagenta+7);

   TLegend* legmc;
   legmc = new TLegend(0.65, 0.7, 0.82, 0.89, "", "brNDC"); 
   legmc->SetTextFont(42);
   legmc->SetBorderSize(0);
   legmc->SetFillStyle(0);
   legmc->AddEntry(h_ntrkJet,"Tracker Jets", "L");
   legmc->AddEntry(h_ngenJet,"Generated Jets", "L");


   h_ntrkJet->Draw("hist"); 
   h_ntrkJet->GetYaxis()->SetTitle("Events");     
   h_ntrkJet->GetXaxis()->SetTitle("Number of Jets");     
   h_ngenJet->Draw("histsame");
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetN.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetN_LOG.png");
  
   //jet pt
   h_trkJetPt->SetLineColor(kAzure+7);
   h_genJetPt->SetLineColor(kMagenta+7);
   h_trkJetPt->Draw("hist"); 
   h_trkJetPt->GetYaxis()->SetTitle("Events");     
   h_trkJetPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");     
   h_genJetPt->Draw("histsame");
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetPt.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetPt_LOG.png");
   
   //jet eta
   h_trkJetEta->SetLineColor(kAzure+7);
   h_genJetEta->SetLineColor(kMagenta+7);
   h_trkJetEta->Draw("hist"); 
   h_trkJetEta->GetYaxis()->SetTitle("Events");     
   h_trkJetEta->GetXaxis()->SetTitle("Jet #eta ");     
   h_genJetEta->Draw("histsame");
   legmc->Draw("same");
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetEta.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetEta_LOG.png");
 
  //jet phi
   h_trkJetPhi->SetLineColor(kAzure+7);
   h_genJetPhi->SetLineColor(kMagenta+7);
   h_trkJetPhi->Draw("hist"); 
   h_trkJetPhi->GetYaxis()->SetTitle("Events");     
   h_trkJetPhi->GetXaxis()->SetTitle("Jet #phi ");     
   h_genJetPhi->Draw("histsame");
   legmc->Draw("same");
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetPhi.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetPhi_LOG.png");
 

   //jet mass
   h_trkJetMass->SetLineColor(kAzure+7);
   h_genJetMass->SetLineColor(kMagenta+7);
   h_trkJetMass->Draw("hist"); 
   h_trkJetMass->GetYaxis()->SetTitle("Events");     
   h_trkJetMass->GetXaxis()->SetTitle("Jet Mass [GeV]");     
   h_genJetMass->Draw("histsame");
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetMass.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetMass_LOG.png");
 

   //jet DR
   h_trkJetDR->SetLineColor(kAzure+7);
   h_trkJetDR->Draw("hist"); 
   h_trkJetDR->GetYaxis()->SetTitle("Events");     
   h_trkJetDR->GetXaxis()->SetTitle("#Delta R (trk-gen)"); 
   c1->SetLogy(0);     
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetDR.png");
   c1->SetLogy();
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/JetDR_LOG.png");
   

   //2D distributions
   c1->SetLogy(0);
   h2_pt_eta->Draw("COLZ");
   h2_pt_eta->GetYaxis()->SetTitle("Tracker Jet p_{T}");     
   h2_pt_eta->GetXaxis()->SetTitle("Tracker Jet #eta");     
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Jet_pt_vs_eta.png");
   
   h2_phi_eta->Draw("COLZ");
   h2_phi_eta->GetYaxis()->SetTitle("Tracker Jet #phi");     
   h2_phi_eta->GetXaxis()->SetTitle("Tracker Jet #eta");     
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Jet_phi_vs_eta.png");
   
   h2_resp_eta->Draw("COLZ");
   h2_resp_eta->GetYaxis()->SetTitle("Tracker Jet Response");     
   h2_resp_eta->GetXaxis()->SetTitle("Tracker Jet #eta");     
   c1->SaveAs("~/www/upgrade/VFPix-trkJets/Jet_response_vs_eta.png");
  
   //means and rms response vs eta
   int etabins=160;
   double mean[etabins];
   double rms[etabins];
   double mean_err[etabins];
   double rms_err[etabins];
 
   for(int k=0; k<10;k+=4){
     std::cout<<k<<std::endl;
     TH1F* h1 = (TH1F*) h2_resp_eta->ProjectionY("",k,k,"e");
     h1->GetXaxis()->SetTitle("<Trk Jet Pt/Gen Jet Pt>");
     h1->GetYaxis()->SetTitle("Entries");
     if(h1->GetEntries()>0) h1->Fit("gaus", "", "",  0., 5.);
     h1->GetFunction("gaus")->SetLineColor(kBlue);
     h1->GetFunction("gaus")->SetLineWidth(2);
     h1->Draw("P");
     mean[k] = h1->GetFunction("gaus")->GetParameter(1);
     mean_err[k] = h1->GetFunction("gaus")->GetParError(1);
     rms[k]  = h1->GetFunction("gaus")->GetParameter(2);
     rms_err[k] = h1->GetFunction("gaus")->GetParError(2);
     c1->SaveAs(TString::Format("~/www/upgrade/VFPix-trkJets/Jet_response_fit_bin_%d.png",k));
   }

}
