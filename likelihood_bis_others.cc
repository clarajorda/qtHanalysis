#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include "THStack.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TLatex.h"
// #include "LatinoPlotOthers.C"
#include "LatinoPlotWithFakes_NewVecbos_all.C"

#include <iostream>
#include <iomanip>

#define NVARIABLES 5
#define NSAMPLES   11
#define LUMINOSITY 19700

#define USECUTJETS 1

TString suffix_="NewVecbos_Nominal_V01_April_cutJets_July";
TString files[NSAMPLES] = {

  "/cmsrm/pc24_2/jorda/data/finalresults_NewVecbos_Nominal_V01_April_cutJets/Data/DataAll.root",

  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/Dibosons-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/Others-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/TTH-Inclusive-125-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/TTWJets-madgraph-v1-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/TTZJets-madgraph-v2-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_v1-datasetAll.root",

  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/qtH-blv_1M-mH125_Ctm1__MLR-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/qtH-blv_1M-mH125_Ctm1__MLR_qtHtautau-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/WtH_mH125_Ct1__MLR_minus1-datasetAll.root",
  "/cmsrm/pc24_2/jorda/data/finalresults_"+suffix_+"/MC/WtH_mH125_Ct1__MLR_minus1_qtHtautau-datasetAll.root"

};

// ---------
// Functions
// ---------

// ! main function
void  likelihood(int mychannel, bool withFakes, bool blindData, bool HistosFromFile, bool systematics, TString modeSystematic, TString tag);

// ! fill the histograms to be used to estimate the likelihood value
void  fillHistogram           (int mychannel, int iSample, bool withFakes, TH1F *hVariable[NVARIABLES]);
void  fillLikelihoodHistogram (int mychannel, int iSample, bool withFakes, bool systematics, TString modeSystematic, TH1F *hVariable_Signal[NVARIABLES], TH1F *hVariable_Background[NVARIABLES], TH1F *hLikelihood);

// ! estimate likelihood for a process
float calculateLikelihood(TH1F *hVariable_Signal[NVARIABLES], TH1F *hVariable_Background[NVARIABLES], double inputs[NVARIABLES]);

// ! set axis titles and colors
void  setHistogram(TH1F *h, TString xtitle, TString ytitle, int color);
void  setHistogram(TH1F *h, TString xtitle, TString ytitle, int color, TString label);

void  fillLikelihoodHistogram (int mychannel, int iSample,  bool withFakes, bool systematics, TString modeSystematic, TH1F *hVariable_Signal[NVARIABLES], TH1F *hVariable_Background[NVARIABLES], TH1F *hLikelihood){

  TString filename = files[iSample];

  std::cout << " :: [fillLikelihoodHistogram Function] :: Sample number [" << iSample << "]. Opening file ... [" << filename << "]" << std::endl;

  TFile *file = TFile::Open(filename);
  file->cd();
  TTree *tree;
  tree = (TTree*)file->Get("latino");
  
  // -------------------------------
  // Tree Variables (reduced format)
  // -------------------------------

  float pfmet;
  float mll;
  float maxPt;
  float pt1;
  float pt2;
  float pt3;
  float ch1,ch2;
  float ht,htjets;
  float minDeltaR;
  int   numCentralMJets;
  int   numTaggedMJets;
  int   numForwardMJets;
  int   numForwardMJets_eta1p7;
  int   numForwardMJets_eta2p0;
  int   numForwardMJets_eta2p2;
  int   numForwardMJets_eta2p4;
  float totalch;
  float dEta_bQuarkMaxTagged_qQuarkMaxEta;
  float dEta_bQuarkMaxPt_qQuarkMaxEta;
  float dEta_bQuarkMaxPt_qQuarkMaxPt;
  float maxEtaForwardJet;

  int isWW;
  int isTauTau;

  float nextra;
  float baseW;
  float puW;
  float effW;
  float btagW,btagWerr;
  
  float channel;

  // fake rate weights
  float        ppf[3],        pff[3],        fff[3];
  float ppfSyst[5][3], pffSyst[5][3], fffSyst[5][3];
  for(int i=0;i<3;i++){
    ppf[i] = 1.0;
    pff[i] = 1.0;
    fff[i] = 1.0;
    for(int j=0;j<5;j++){
      ppfSyst[j][i] = 1.0;
      pffSyst[j][i] = 1.0;
      fffSyst[j][i] = 1.0;
    }
  }  

  // ------------------------
  // Assign variables in tree
  // ------------------------
 
  tree->SetBranchAddress("pfmet", &pfmet);
  tree->SetBranchAddress("mll", &mll);
  tree->SetBranchAddress("maxPt", &maxPt);
  tree->SetBranchAddress("pt1", &pt1);
  tree->SetBranchAddress("pt2", &pt2);
  tree->SetBranchAddress("pt3", &pt3);
  tree->SetBranchAddress("ch1", &ch1);
  tree->SetBranchAddress("ch2", &ch2);
  tree->SetBranchAddress("ht", &ht);
  tree->SetBranchAddress("htjets", &htjets);
  tree->SetBranchAddress("minDeltaR", &minDeltaR);
  tree->SetBranchAddress("numCentralMJets", &numCentralMJets);
  tree->SetBranchAddress("numTaggedMJets", &numTaggedMJets);
  tree->SetBranchAddress("numForwardMJets", &numForwardMJets);
  tree->SetBranchAddress("numForwardMJets_eta1p7", &numForwardMJets_eta1p7);
  tree->SetBranchAddress("numForwardMJets_eta2p0", &numForwardMJets_eta2p0);
  tree->SetBranchAddress("numForwardMJets_eta2p2", &numForwardMJets_eta2p2);
  tree->SetBranchAddress("numForwardMJets_eta2p4", &numForwardMJets_eta2p4);
  tree->SetBranchAddress("totalch", &totalch);
  tree->SetBranchAddress("dEta_bQuarkMaxTagged_qQuarkMaxEta", &dEta_bQuarkMaxTagged_qQuarkMaxEta);
  tree->SetBranchAddress("dEta_bQuarkMaxPt_qQuarkMaxEta", &dEta_bQuarkMaxPt_qQuarkMaxEta);
  tree->SetBranchAddress("dEta_bQuarkMaxPt_qQuarkMaxPt",&dEta_bQuarkMaxPt_qQuarkMaxPt);
  tree->SetBranchAddress("maxEtaForwardJet", &maxEtaForwardJet);

  if(iSample==itHWW || iSample==itHtt || iSample==iWtHWW || iSample==iWtHtt){
    tree->SetBranchAddress("isWW",&isWW);
    tree->SetBranchAddress("isTauTau",&isTauTau);
    tree->SetBranchAddress("btagW", &btagW);
    tree->SetBranchAddress("btagWerr", &btagWerr);
  }

  tree->SetBranchAddress("nextra", &nextra);
  tree->SetBranchAddress("baseW", &baseW);
  tree->SetBranchAddress("puW", &puW);
  tree->SetBranchAddress("effW", &effW);

  tree->SetBranchAddress("channel", &channel);

  if (withFakes && iSample==itt){
    tree->SetBranchAddress("ppf",&ppf);
    tree->SetBranchAddress("pff",&pff);
    tree->SetBranchAddress("fff",&fff);
    tree->SetBranchAddress("ppfSyst",&ppfSyst);
    tree->SetBranchAddress("pffSyst",&pffSyst);
    tree->SetBranchAddress("fffSyst",&fffSyst);
  }

  // ----
  // Loop
  // ----

  file->cd();
  int n = tree->GetEntries();

  double inputs[NVARIABLES];
  Float_t likelihood = -1;

  // Create also a new ROOTFILE with the Likelihood variable
  TFile *newFile = new TFile(filename+"_likelihood","recreate");
  TTree *newTree = new TTree("latino","latino tree with likelihood variable");


  newTree->Branch("baseW", &baseW, "baseW/F");
  newTree->Branch("puW"  , &puW  , "puW/F");
  newTree->Branch("effW" , &effW , "effW/F");

  newTree->Branch("btagW"   , &btagW   , "btagW/F");
  newTree->Branch("btagWerr", &btagWerr, "btagWerr/F");

  newTree->Branch("channel", &channel, "channel/F");
  newTree->Branch("nextra", &nextra, "nextra/F");

  newTree->Branch("pfmet", &pfmet, "pfmet/F");
  newTree->Branch("mll", &mll, "mll/F");
  newTree->Branch("pt1", &pt1, "pt1/F");
  newTree->Branch("pt2", &pt2, "pt2/F");
  newTree->Branch("pt3", &pt3, "pt3/F");
  newTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/F");
  newTree->Branch("numTaggedMJets", &numTaggedMJets, "numTaggedMJets/I");
  newTree->Branch("numForwardMJets", &numForwardMJets, "numForwardMJets/I");
  newTree->Branch("totalch", &totalch, "totalch/F");
  newTree->Branch("dEta_bQuarkMaxTagged_qQuarkMaxEta", &dEta_bQuarkMaxTagged_qQuarkMaxEta, "dEta_bQuarkMaxTagged_qQuarkMaxEta/F");
  newTree->Branch("likelihood", &likelihood, "likelihood/F");

  if (systematics && (modeSystematic=="Fakes_met" || modeSystematic=="Fakes_dphiL" || modeSystematic=="Fakes_dphiH" || modeSystematic=="Fakes_ptTagL" || modeSystematic=="Fakes_ptTagH") ) 
    std::cout << "[SYSTEMATICS SET] ... on FAKES, mode: " << modeSystematic << std::endl;

  if (systematics && (modeSystematic=="BTag_Up" || modeSystematic=="BTag_Down"))
    std::cout << "[SYSTEMATICS SET] ... on B-TAGGING, mode: " << modeSystematic << std::endl;

  // event loop
  for(int e=0;e<n;e++){
  // for(int e=0;e<500;e++){
    
    tree->GetEntry(e);
    
    // cuts for the event selection
    if ( nextra == 0 && maxPt > 20. && pfmet > 30. && mll > 20. && numForwardMJets > 0 && numTaggedMJets == 1 && ( ((channel==2 || channel==3) && ch1*ch2>0) || (mll<75. || mll>105.) ) ){ // LIKELIHOOD 

      // if ( nextra == 0 && maxPt > 20. && pfmet > 30. && mll > 20. && numForwardMJets < 1 && numTaggedMJets == 1 && ( ((channel==2 || channel==3) && ch1*ch2>0) || (mll<75. || mll>105.) ) ){ // LIKELIHOOD 


      if( ( mychannel!=channel && mychannel!=-1) ) continue;

      // if (numForwardMJets_eta2p4 > 0) continue; // NOW

      // weight for the event, taking into account the cross-section
      float weight = 1.0;
      if ( iSample != iData && iSample != itt && iSample != itHWW && iSample != iWtHWW && iSample != itHtt && iSample != iWtHtt ){
	weight = baseW*puW*effW*LUMINOSITY;

      }else if(iSample == itt){
	if(!withFakes){ // MC
	  weight = baseW*puW*effW*LUMINOSITY;
	}else if (withFakes){ // Data-driven
	  if      (!systematics){
	    weight = weight*(ppf[0] + pff[0] + fff[0]);

	  }else if(systematics){

	    if (modeSystematic=="Fakes_StatFR_Up"){
	      weight = weight*(ppf[1] + pff[1] + fff[1]);
	    }
	    else if (modeSystematic=="Fakes_StatFR_Down"){
	      weight = weight*(ppf[2] + pff[2] + fff[2]);
	    }
	    else if (modeSystematic=="Fakes_met"){
	      weight = weight*(ppfSyst[0][0] + pffSyst[0][0] + fffSyst[0][0]);
	    }
	    else if (modeSystematic=="Fakes_dphiL"){
	      weight = weight*(ppfSyst[1][0] + pffSyst[1][0] + fffSyst[1][0]);
	    }
	    else if (modeSystematic=="Fakes_dphiH"){
	      weight = weight*(ppfSyst[2][0] + pffSyst[2][0] + fffSyst[2][0]);
	    }
	    else if (modeSystematic=="Fakes_ptTagL"){
	      weight = weight*(ppfSyst[3][0] + pffSyst[3][0] + fffSyst[3][0]);
	    }
	    else if (modeSystematic=="Fakes_ptTagH"){
	      weight = weight*(ppfSyst[4][0] + pffSyst[4][0] + fffSyst[4][0]);
	    }
	    else{
	      weight = weight*(ppf[0] + pff[0] + fff[0]);
	    }

	  }
  
	}
      }

      if      (iSample == itHWW || iSample == iWtHWW ) {
	if   (isWW==1) {  weight = baseW*puW*effW*LUMINOSITY*isWW; }
	else           {  continue; }
      }

      else if (iSample == itHtt || iSample == iWtHtt) {
	if   (isTauTau==1) { weight = baseW*puW*effW*LUMINOSITY*isTauTau; }
	else               {  continue; }
      }

      if ( (iSample == itHWW || iSample == itHtt || iSample == iWtHWW || iSample == iWtHtt) ){

	if (!systematics){
	  weight = weight * btagW;	  
	}else{
	  if      (modeSystematic=="BTag_Up"  ) { weight = weight * (btagW+btagWerr); }
	  else if (modeSystematic=="BTag_Down") { weight = weight * (btagW-btagWerr); }
	  else                                  { weight = weight * btagW; }
	}
	
      }

      // give variables used for the likelihood training, *in the same order*
      // inputs[0] = numTaggedMJets;
      inputs[0] = numCentralMJets;
      inputs[1] = numForwardMJets_eta2p4;
      inputs[2] = minDeltaR;
      inputs[3] = totalch;

      // inputs[4] = abs(dEta_bQuarkMaxPt_qQuarkMaxPt); // A
      inputs[4] = abs(dEta_bQuarkMaxPt_qQuarkMaxEta); // B
      // inputs[4] = abs(maxEtaForwardJet); // C

      // inputs[5] = mll;
      // inputs[5] = htjets;

      // inputs[6] = pfmet;

      likelihood = calculateLikelihood(hVariable_Signal,hVariable_Background,inputs);
      hLikelihood -> Fill(likelihood, weight);

      newTree->Fill();      
      
    }// loop events
    
  }

  newFile->cd();
  newTree->Write();      
  newFile->Close();

  file->cd();
  file->Close();

}

void  fillHistogram(int mychannel, int iSample, bool withFakes, TH1F *hVariable[NVARIABLES]){

  TString filename = files[iSample];

  TFile *file = TFile::Open(filename);

  std::cout << " :: [fillHistogram Function] :: Sample number [" << iSample << "]. Opening file ... [" << filename << "]" << std::endl;

  TTree *tree;
  tree = (TTree*)file->Get("latino");

  // -------------------------------
  // Tree Variables (reduced format)
  // -------------------------------
  float pfmet;
  float mll;
  float maxPt;
  float pt1;
  float pt2;
  float pt3;
  float ch1,ch2;
  float ht,htjets;
  float minDeltaR;
  int   numCentralMJets;
  int   numTaggedMJets;
  int   numForwardMJets;
  int   numForwardMJets_eta1p7;
  int   numForwardMJets_eta2p0;
  int   numForwardMJets_eta2p2;
  int   numForwardMJets_eta2p4;
  float totalch;
  float dEta_bQuarkMaxTagged_qQuarkMaxEta;
  float dEta_bQuarkMaxPt_qQuarkMaxEta;
  float dEta_bQuarkMaxPt_qQuarkMaxPt;
  float maxEtaForwardJet;

  int isWW;
  int isTauTau;

  float nextra;
  float baseW;
  float puW;
  float effW;
  float btagW;

  float channel;

  // fake rate weights
  float ppf[3], pff[3], fff[3];
  for(int i=0;i<3;i++){
    ppf[i] = 1.0;
    pff[i] = 1.0;
    fff[i] = 1.0;
  }


  // ------------------------
  // Assign variables in tree
  // ------------------------
 
  tree->SetBranchAddress("pfmet", &pfmet);
  tree->SetBranchAddress("mll", &mll);
  tree->SetBranchAddress("maxPt", &maxPt);
  tree->SetBranchAddress("pt1", &pt1);
  tree->SetBranchAddress("pt2", &pt2);
  tree->SetBranchAddress("pt3", &pt3);
  tree->SetBranchAddress("ch1", &ch1);
  tree->SetBranchAddress("ch2", &ch2);
  tree->SetBranchAddress("ht", &ht);
  tree->SetBranchAddress("htjets", &htjets);
  tree->SetBranchAddress("minDeltaR", &minDeltaR);
  tree->SetBranchAddress("numCentralMJets", &numCentralMJets);
  tree->SetBranchAddress("numTaggedMJets", &numTaggedMJets);
  tree->SetBranchAddress("numForwardMJets", &numForwardMJets);
  tree->SetBranchAddress("numForwardMJets_eta1p7", &numForwardMJets_eta1p7);
  tree->SetBranchAddress("numForwardMJets_eta2p0", &numForwardMJets_eta2p0);
  tree->SetBranchAddress("numForwardMJets_eta2p2", &numForwardMJets_eta2p2);
  tree->SetBranchAddress("numForwardMJets_eta2p4", &numForwardMJets_eta2p4);
  tree->SetBranchAddress("totalch", &totalch);
  tree->SetBranchAddress("dEta_bQuarkMaxTagged_qQuarkMaxEta", &dEta_bQuarkMaxTagged_qQuarkMaxEta);
  tree->SetBranchAddress("dEta_bQuarkMaxPt_qQuarkMaxEta", &dEta_bQuarkMaxPt_qQuarkMaxEta);
  tree->SetBranchAddress("dEta_bQuarkMaxPt_qQuarkMaxPt",&dEta_bQuarkMaxPt_qQuarkMaxPt);
  tree->SetBranchAddress("maxEtaForwardJet", &maxEtaForwardJet);
  
  if(iSample==itHWW || iSample==itHtt || iSample==iWtHWW || iSample==iWtHtt){
    tree->SetBranchAddress("isWW",&isWW);
    tree->SetBranchAddress("isTauTau",&isTauTau);
    tree->SetBranchAddress("btagW",&btagW);
  }

  tree->SetBranchAddress("nextra", &nextra);
  tree->SetBranchAddress("baseW", &baseW);
  tree->SetBranchAddress("puW", &puW);
  tree->SetBranchAddress("effW", &effW);

  tree->SetBranchAddress("channel", &channel);

  if(withFakes && iSample==itt){
    tree->SetBranchAddress("ppf",&ppf);
    tree->SetBranchAddress("pff",&pff);
    tree->SetBranchAddress("fff",&fff);
  }

  // ----
  // Loop
  // ----

  int n = tree->GetEntries();

  // event loop
  for(int e=0;e<n;e++){
    
    tree->GetEntry(e);

    // Likelihood
    if ( nextra == 0 && maxPt > 20. && pfmet > 30. && mll > 20. && numForwardMJets > 0 && numTaggedMJets == 1 && ( ((channel==2 || channel==3) && ch1*ch2>0) || (mll<75. || mll>105.) ) ){ // LIKELIHOOD 

      // if ( nextra == 0 && maxPt > 20. && pfmet > 30. && mll > 20. && numForwardMJets < 1 && numTaggedMJets == 1 && ( ((channel==2 || channel==3) && ch1*ch2>0) || (mll<75. || mll>105.) ) ){ // LIKELIHOOD 


      if( ( mychannel!=channel && mychannel!=-1 ) ) continue;

      // weight for the event, taking into account the cross-section
      float weight = 1.0;

      if ( iSample != iData && iSample != itt && iSample != itHWW && iSample != iWtHWW && iSample != itHtt && iSample != iWtHtt){
	weight = baseW*puW*effW*LUMINOSITY;

      }else if(iSample == itt){
	if(!withFakes) // MC
	  weight = baseW*puW*effW*LUMINOSITY;

	else if (withFakes) //Data-driven
	  weight = weight*(ppf[0] + pff[0] + fff[0]);
      }

      if      (iSample == itHWW || iSample == iWtHWW) {
	if   (isWW==1) { weight = baseW*puW*effW*LUMINOSITY * isWW * btagW; }
	else           { continue; }
      }

      else if (iSample == itHtt || iSample == iWtHtt) {
	if   (isTauTau==1) { weight = baseW*puW*effW*LUMINOSITY * isTauTau * btagW; }
	else               { continue; }
      }

      hVariable[0] -> Fill( numCentralMJets, weight);
      hVariable[1] -> Fill( numForwardMJets_eta2p4, weight);
      hVariable[2] -> Fill( minDeltaR, weight );
      hVariable[3] -> Fill( totalch, weight);

      // hVariable[4] -> Fill( abs(dEta_bQuarkMaxPt_qQuarkMaxPt), weight); // A      
      hVariable[4] -> Fill( abs(dEta_bQuarkMaxPt_qQuarkMaxEta), weight); // B
      // hVariable[4] -> Fill( abs(maxEtaForwardJet), weight ); // C

      // hVariable[5] -> Fill(   mll, weight );
      // hVariable[5] -> Fill(   htjets, weight );

      // hVariable[6] -> Fill( pfmet, weight );

    }

  } // loop events

}// fillHistogram

// Define Functions
void setHistogram(TH1F *h, TString xtitle, TString ytitle, int color, TString label){

   // h->SetNormFactor(1.0);

  h->SetTitle(label);
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(color);
  h->SetLineWidth(2);

}

void setHistogram(TH1F *h, TString xtitle, TString ytitle, int color){

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->SetLineColor(color);
  h->SetLineWidth(2);

}

// ! main function
void likelihood(int mychannel, bool withFakes, bool blindData, bool HistosFromFile, bool systematics, TString modeSystematic, TString tag){

  std::cout << "estoy aqui" << std::endl;

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x LatinoStyle.C");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetMarkerColor(1);

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  if (withFakes){
    if(USECUTJETS)
      files[itt]="/cmsrm/pc24_2/jorda/data/finalresults_NewVecbos_FakeRateApplication_V01_Data_April_cutJets/Data/DataAll.root";
    else
      files[itt]="/cmsrm/pc24_2/jorda/data/finalresults_NewVecbos_FakeRateApplication_V01_Data_April_mvaJets/Data/DataAll.root";
  }
  // ---------------------------
  // Define processes to be used
  // ---------------------------

  TString species[NSAMPLES];
  species[0]="Data";
  species[1]="VV";
  species[2]="Others";
  species[3]="ttH";
  species[4]="ttW";
  species[5]="ttZ";
  species[6]="Fakes";
  species[7]="tHWW";
  species[8]="tHtt"; 
  species[9]="WtHWW";
  species[10]="WtHtt"; 

  TString variables[NVARIABLES];
  // variables[0] = "numTaggedMJets";
  variables[0] = "numCentralMJets";
  // variables[1] = "numForwardMJets"; // no
  variables[1] = "numForwardMJets_eta2p4";
  variables[2] = "minDeltaR";
  variables[3] = "totalch";
  // variables[4] = "dEta_bQuarkMaxPt_qQuarkMaxPt"; // A
  variables[4] = "dEta_bQuarkMaxPt_qQuarkMaxEta"; // B
  // variables[4] = "maxEtaForwardJet"; // C
  // variables[5] = "mll";
  // variables[5] = "htjets";
  // variables[6] = "pfmet";

  TString xaxisLabel[NVARIABLES];
  // xaxisLabel[0] = "#b-tagged CVS-M jets";
  xaxisLabel[0] = "#central jets";
  xaxisLabel[1] = "#forward jets |\\eta| > 2.4";
  xaxisLabel[2] = "min (\\Delta R(l,l))";
  xaxisLabel[3] = "total charge";
  xaxisLabel[4] = "\\Delta\\eta(q,b)";
  // xaxisLabel[4] = "max(\\eta^{forward jet})";
  // xaxisLabel[5] = "m_{ll}";
  // xaxisLabel[5] = "H_{T} jets";
  // xaxisLabel[6] = "PF E_{T}^{miss}";
  
  TString units[NVARIABLES];
  units[0] = "";
  units[1] = "";
  units[2] = "";
  units[3] = "";
  units[4] = "";
  // units[5] = "GeV";
  // units[6] = "GeV";

  int nbins[NVARIABLES];
  nbins[0] = 8;  
  nbins[1] = 8;  
  nbins[2] = 15;
  nbins[3] = 6;  
  nbins[4] = 12;

  // nbins[5] = 18; //20
  // nbins[5] = 20; //20
  // nbins[6] = 14; //15

  float range[NVARIABLES][2]; // N variables, min, max
  // njets
  range[0][0] = 0.;
  range[0][1] = 8.; // 8. 
  // njets
  range[1][0] = 0.;
  range[1][1] = 8.; // 8.
  // mindR
  range[2][0] = 0.;
  range[2][1] = 5.; // 5.  
  // total q
  range[3][0] = -2;
  range[3][1] = +2;
  // delta eta
  //range[4][0] =  0.0;
  //range[4][1] = +7.0; // 7.

  // max forward eta
  range[4][0] = 0.;
  range[4][1] = 5.; // 5.

  // mll
  //range[5][0] = 20.;
  //range[5][1] = 200.; // 1000.

  //range[5][0] = 20.;
  //range[5][1] = 600.; // 1000.

  //met
  // range[6][0] = 30.;
  // range[6][1] = 250.;

  // ----------
  // Histograms
  // ----------

  // variables
  TH1F *hVariables[NSAMPLES][NVARIABLES]; 

  if (HistosFromFile){

    std::cout << "[Individual Histograms] I am READING the histograms from a file " << std::endl;

    TFile* histosFile;

    if (withFakes)
      histosFile = TFile::Open("/afs/cern.ch/work/j/jorda/prueba4/CMSSW_5_3_6/src/TopHiggsAnalysis/Likelihood/LikelihoodHISTOS/Likelihood_histos_withFakes_PreviousHistograms_"+tag+".root");
    else
      histosFile = TFile::Open("/afs/cern.ch/work/j/jorda/prueba4/CMSSW_5_3_6/src/TopHiggsAnalysis/Likelihood/LikelihoodHISTOS/Likelihood_histos_withMC_PreviousHistograms_"+tag+".root");

    for (int s=0;s<NSAMPLES;s++){
      for (int v=0;v<NVARIABLES;v++){

	hVariables[s][v] = (TH1F*) histosFile->Get("hVariables_"+variables[v]+"_"+species[s])->Clone("hVariables_"+variables[v]+"_"+species[s]);
	hVariables[s][v]->SetDirectory(0);

	if (v == 0) std::cout << "Sample [" << s << "]. Histogram entries [" << hVariables[s][v]->GetEntries() <<"]" << std::endl;

      }
    }

  }else{

    std::cout << "[Individual Histograms] I am PRODUCING the histograms " << std::endl;

    for (int s=0;s<NSAMPLES;s++){
      for (int v=0;v<NVARIABLES;v++){
	hVariables[s][v] = new TH1F("hVariables_"+variables[v]+"_"+species[s], "", nbins[v], range[v][0], range[v][1] );
	hVariables[s][v] -> Sumw2();
      }
    }
    
    // ! fill the histograms
    for (int s=0;s<NSAMPLES;s++){
      fillHistogram (mychannel, s, withFakes, hVariables[s] ); 
    }

  }

  TH1F *hVariablesTotalSignal    [NVARIABLES];   
  TH1F *hVariablesTotalBackground[NVARIABLES];   

  if (HistosFromFile){

    std::cout << "[Total S and B] I am READING the histograms from a file " << std::endl;

    TFile* histosFile;
    
    if (withFakes)
      histosFile = TFile::Open("/afs/cern.ch/work/j/jorda/prueba4/CMSSW_5_3_6/src/TopHiggsAnalysis/Likelihood/LikelihoodHISTOS/Likelihood_histos_withFakes_PreviousHistograms_"+tag+".root");
    else
      histosFile = TFile::Open("/afs/cern.ch/work/j/jorda/prueba4/CMSSW_5_3_6/src/TopHiggsAnalysis/Likelihood/LikelihoodHISTOS/Likelihood_histos_withMC_PreviousHistograms_"+tag+".root");
    
    for (int v=0;v<NVARIABLES;v++){
      
      hVariablesTotalSignal[v] = (TH1F*) histosFile->Get("hVariablesTotalSignal_"+variables[v])->Clone("hVariablesTotalSignal_"+variables[v]);
      hVariablesTotalSignal[v]->SetDirectory(0);
      
      hVariablesTotalBackground[v] = (TH1F*) histosFile->Get("hVariablesTotalBackground_"+variables[v])->Clone("hVariablesTotalBackground_"+variables[v]);
      hVariablesTotalBackground[v]->SetDirectory(0);
            
    }
        
  }else{

    std::cout << "[Total S and B] I am PRODUCING the histograms " << std::endl;

    for (int v=0;v<NVARIABLES;v++){
      
      // total signal histogram
      hVariablesTotalSignal[v] = new TH1F("hVariablesTotalSignal_"+variables[v], "", nbins[v], range[v][0], range[v][1] );
      hVariablesTotalSignal[v] -> Sumw2();
      
      // total signal tH(WW) and tH(tautau)
      hVariablesTotalSignal[v] -> Add(hVariables[itHWW][v]);
      hVariablesTotalSignal[v] -> Add(hVariables[itHtt][v]);
      
      hVariablesTotalSignal[v] -> Scale( 1./(hVariablesTotalSignal[v]->Integral()) );
      
      // total background histogram
      hVariablesTotalBackground[v] = new TH1F("hVariablesTotalBackground_"+variables[v], "", nbins[v], range[v][0], range[v][1] );
      hVariablesTotalBackground[v] -> Sumw2();
      
      // total background is WZ + tt + ttw + ttz
      // hVariablesTotalBackground[v] -> Add(hVariables[iVV][v]);
      hVariablesTotalBackground[v] -> Add(hVariables[itt][v]);
      hVariablesTotalBackground[v] -> Add(hVariables[ittw][v]);
      hVariablesTotalBackground[v] -> Add(hVariables[ittz][v]);
      
      hVariablesTotalBackground[v] -> Scale( 1./(hVariablesTotalBackground[v]->Integral()) );
      
    }

  }

  
  // ! normalize to unity
  for (int s=0;s<NSAMPLES;s++){
    for (int v=0;v<NVARIABLES;v++){
      hVariables[s][v] -> Scale( 1./(hVariables[s][v]->Integral()) );
    }
  } 
  // estimate likelihood


  int nBinsLikelihood = 10;
  TH1F *hLikelihood[NSAMPLES];
  TH1F *hLikelihoodBlinded;

  if (systematics){
    
    if (!HistosFromFile){
      std::cout << "[SYSTEMATICS SET] but PRODUCING histos. Exiting!! " << std::endl;
      return;
    }
    
    std::cout << "[SYSTEMATICS SET] Taking flag :: " << modeSystematic << std::endl;
    
    
    
    // For all
    if ( modeSystematic=="MuonMomentum_Yes" || modeSystematic=="ElectronMomentum_Up" || modeSystematic=="ElectronMomentum_Down" || modeSystematic=="JESUncertainty_Up" || 
	 modeSystematic=="JESUncertainty_Down" || modeSystematic=="LeptonSF_Up" || modeSystematic=="LeptonSF_Down" || modeSystematic=="LeptonSF_Yes"){

      // Loop over all samples
      for(int f=1; f<NSAMPLES;f++){
	
	if (f==itt) continue;
	
	files[f].ReplaceAll("Nominal_V01_April","Systematic_"+modeSystematic+"_V01_May");
	std::cout << "[SYSTEMATICS SET] Using for sample : " << f << ", filename : " << files[f] << std::endl;
	
      }
      
    }else if( modeSystematic=="Fakes_met" || modeSystematic=="Fakes_dphiL" || modeSystematic=="Fakes_dphiH" || modeSystematic=="Fakes_ptTagL" || modeSystematic=="Fakes_ptTagH" || 
	      modeSystematic=="Fakes_StatFR_Up" || modeSystematic=="Fakes_StatFR_Down" ){
      std::cout << "[SYSTEMATICS SET] ... on FAKES " << std::endl;

    }else if( modeSystematic=="BTag_Up" || modeSystematic=="BTag_Down" ){
      std::cout << "[SYSTEMATICS SET] ... on B-TAGGING " << std::endl;

    }else{
      std::cout << "[SYSTEMATICS SET] WARNING !! You didn't give any mode! Exiting!! " << std::endl;
      return;
    }   

  }
  

  hLikelihoodBlinded = new TH1F("hLikelihood_Data_Blinded","", (nBinsLikelihood/2), 0.0, 0.5);
  fillLikelihoodHistogram (mychannel, 0, withFakes, systematics, modeSystematic, hVariables[itHWW], hVariablesTotalBackground, hLikelihoodBlinded);

  for (int s=0;s<NSAMPLES;s++){
    
    if (blindData && s == iData){
      hLikelihood[s] = new TH1F("hLikelihood_"+species[s], "", nBinsLikelihood, 0.0, 1.0);
      // hLikelihood[s] = new TH1F("hLikelihood_"+species[s], "", (nBinsLikelihood/2), 0.0, 0.5);
      // hLikelihood[s] = new TH1F("hLikelihood_"+species[s], "", 10, 0.0, 0.5);
    }else{
      hLikelihood[s] = new TH1F("hLikelihood_"+species[s], "", nBinsLikelihood, 0.0, 1.0);
    }

    hLikelihood[s] -> Sumw2();
    
    // CJ nominal
    // fillLikelihoodHistogram (mychannel, s, withFakes, systematics, modeSystematic, hVariablesTotalSignal, hVariablesTotalBackground, hLikelihood[s]);
    fillLikelihoodHistogram (mychannel, s, withFakes, systematics, modeSystematic, hVariables[itHWW], hVariablesTotalBackground, hLikelihood[s]);
    
  }

  LatinoPlot* myPlot = new LatinoPlot();
  myPlot->setLumi(LUMINOSITY/1000.);
  myPlot->addLabel("");
  myPlot->setLabel("likelihood");
  myPlot->setUnits("");

  myPlot->setMCHist(iVV,  hLikelihood[iVV]);
  myPlot->setMCHist(iOthers, hLikelihood[iOthers]);
  myPlot->setMCHist(itth,  hLikelihood[itth]);
  myPlot->setMCHist(ittw, hLikelihood[ittw]);
  myPlot->setMCHist(ittz, hLikelihood[ittz]);
  myPlot->setMCHist(itt,  hLikelihood[itt]);
  myPlot->setMCHist(itHWW, hLikelihood[itHWW]);
  myPlot->setMCHist(itHtt, hLikelihood[itHtt]);
  myPlot->setMCHist(iWtHWW, hLikelihood[iWtHWW]);
  myPlot->setMCHist(iWtHtt, hLikelihood[iWtHtt]);

  if(!blindData) 
    myPlot->setDataHist(hLikelihood[iData]);
  else
    myPlot->setDataHist(hLikelihoodBlinded);

  TString suffix = "";
  if (withFakes) suffix = "withFakes";
  else suffix = "withMC";

  // create the folder
  TString createFolder = "mkdir plots_Likelihood_kk_"+tag+" > /dev/null";
  // system ((char*) createFolder);
  system (createFolder);
  
  TCanvas *c1 = new TCanvas("Likelihood_lin","Likelihood_lin");
  c1->cd();
  myPlot->Draw(1);
  c1->GetFrame()->DrawClone();
  c1->SaveAs("plots_Likelihood_kk_"+tag+"/Likelihood_"+suffix+".lin.png");
  c1->SaveAs("plots_Likelihood_kk_"+tag+"/Likelihood_"+suffix+"_lin.pdf");

  TCanvas *c2 = new TCanvas("Likelihood_log","Likelihood_log");
  c2->cd();
  c2->SetLogy(1);
  myPlot->Draw(1);
  c2->GetFrame()->DrawClone();
  c2->SaveAs("plots_Likelihood_kk_"+tag+"/Likelihood_"+suffix+".log.png");
  c2->SaveAs("plots_Likelihood_kk_"+tag+"/Likelihood_"+suffix+"_log.pdf");

  LatinoPlot* myPlotVariables[NVARIABLES] ;
  TCanvas* canvas[NVARIABLES];

  for(int p=0;p<NVARIABLES;p++){

    canvas[p] = new TCanvas(variables[p],variables[p]);

    myPlotVariables[p] = new LatinoPlot();
    myPlotVariables[p]->setLumi(19.7);

    myPlotVariables[p]->setLabel((xaxisLabel[p]).Data());

    // myPlotVariables[p]->setMCHist(itHWW, hVariables[itHWW][p]);
    myPlotVariables[p]->setMCHist(itHWW, hVariablesTotalSignal[p]);

    // myPlotVariables[p]->setMCHist(itt  , hVariables  [itt][p]);
    myPlotVariables[p]->setMCHist(itt  , hVariablesTotalBackground[p]);

    canvas[p]->cd();
    canvas[p]->SetLogy(1); // Log
    myPlotVariables[p]->Draw(1);
    canvas[p]->GetFrame()->DrawClone();
    //canvas[p]->SaveAs("plots_Likelihood_kk_"+tag+"/PDF_"+variables[p]+"_"+suffix+"_lin.png");
    //canvas[p]->SaveAs("plots_Likelihood_kk_"+tag+"/PDF_"+variables[p]+"_"+suffix+"_lin.pdf");
    canvas[p]->SaveAs("plots_Likelihood_kk_"+tag+"/PDF_"+variables[p]+"_"+suffix+"_log.png");
    canvas[p]->SaveAs("plots_Likelihood_kk_"+tag+"/PDF_"+variables[p]+"_"+suffix+"_log.pdf");
    delete myPlotVariables[p];

  }

  TCanvas *auxCanvas = new TCanvas("auxCanvas","auxCanvas");
  auxCanvas->cd();
  // Draw plots with a dummy LatinoPlot, to set histograms
  LatinoPlot* dummy[NVARIABLES] ;
  for(int p=0;p<NVARIABLES;p++){

    dummy[p] = new LatinoPlot();
    dummy[p]->setLumi(19.7);

    dummy[p]->setLabel((xaxisLabel[p]).Data());

    dummy[p]->setMCHist(itHWW, hVariables[itHWW][p]);
    dummy[p]->setMCHist(itHtt, hVariables[itHtt][p]);
    dummy[p]->setMCHist(iWtHWW, hVariables[iWtHWW][p]);
    dummy[p]->setMCHist(iWtHtt, hVariables[iWtHtt][p]);
    dummy[p]->setMCHist( itt, hVariables[ itt][p]);
    dummy[p]->setMCHist( iVV, hVariables[ iVV][p]);
    dummy[p]->setMCHist(ittw, hVariables[ittw][p]);
    dummy[p]->setMCHist(ittz, hVariables[ittz][p]);
    dummy[p]->setMCHist(iOthers, hVariables[iOthers][p]);
    dummy[p]->setMCHist(itth, hVariables[itth][p]);

    dummy[p]->Draw(1);
    delete dummy[p];

  }

  float cut = 0.;
  int binCut = hLikelihood[itHWW]->FindBin(cut);
  std::cout << " Cut for Likelihood is : " << cut << " corresponding to the bin " << binCut << " in the histogram " << std::endl;

  double nEvents [NSAMPLES]; // expected events
  double neEvents[NSAMPLES]; // error expected events

  std::cout<<"  "<< std::setiosflags(std::ios::fixed) << std::setprecision(3) << endl;

  for (int s=0; s<NSAMPLES; s++){

    // expected number of events and error
    nEvents [s] = hLikelihood[s]->Integral();
    double entries = hLikelihood[s]->GetEntries();
    neEvents[s] = (entries>0)  ? (nEvents[s]*TMath::Sqrt(entries)/entries) : 0.0;

    if (blindData && s==iData) continue;

    std::cout << " " << species[s] << " = (" << nEvents [s] << " +- " << neEvents[s] << ")" << std::endl; 

  }

  float expbkg = nEvents[iVV] + nEvents[itt] + nEvents[ittw] + nEvents[ittz] + nEvents[iOthers] + nEvents[itth];
  float errbkg = TMath::Sqrt( neEvents[iVV]*neEvents[iVV] + neEvents[itt]*neEvents[itt] + neEvents[ittw]*neEvents[ittw] + neEvents[ittz]*neEvents[ittz] + 
			      neEvents[iOthers]*neEvents[iOthers] + neEvents[itth]*neEvents[itth]);

  std::cout << " " << std::endl;
  std::cout << " Total background = ( " << expbkg << " +- " << errbkg << " ) " << std::endl;
  std::cout<<"  "<< std::setiosflags(std::ios::fixed) << std::setprecision(6) << endl; 
  
  TFile *fOut;
  if (withFakes){
    if(!systematics)
      fOut = new TFile("LikelihoodHISTOS/Likelihood_histos_withFakes_PreviousHistograms_"+tag+".root","RECREATE");
    else
      fOut = new TFile("LikelihoodHISTOS/Likelihood_histos_withFakes_PreviousHistograms_Systematics_"+modeSystematic+"_"+tag+".root","RECREATE");
  }  else {
    fOut = new TFile("LikelihoodHISTOS/Likelihood_histos_withMC_PreviousHistograms_"+tag+".root","RECREATE");
  }
  fOut->cd();

  for (int s=0;s<NSAMPLES;s++){
    hLikelihood[s]->Write();  
    for(int p=0;p<NVARIABLES;p++){
      hVariables[s][p]->Write();  
    }
  }

  hLikelihoodBlinded->Write();

  for(int p=0;p<NVARIABLES;p++){
    hVariablesTotalSignal    [p]->Write();
    hVariablesTotalBackground[p]->Write();
  }

  fOut->Write();
  fOut->Close();

}// void likelihood()

// ! estimate likelihood for a process
float calculateLikelihood(TH1F *hVariable_Signal[NVARIABLES], TH1F *hVariable_Background[NVARIABLES], double inputs[NVARIABLES]){

  float likelihood_s = 1.;
  float likelihood_b = 1.;
  
  for(int i=0; i<NVARIABLES; i++){
  
    // if the variable lays in the overflow bin, then take last bin from the histogram
    int last_bin_s = hVariable_Signal    [i]->FindLastBinAbove(-100.,1);
    int last_bin_b = hVariable_Background[i]->FindLastBinAbove(-100.,1);

    int bin_s = hVariable_Signal[i]     -> FindBin(inputs[i]);
    int bin_b = hVariable_Background[i] -> FindBin(inputs[i]);

    /*
    if (bin_s > last_bin_s){
      std::cout << " i variable is " << i << " , value input " << inputs[i] << std::endl;
      std::cout << " bin_s is higher than last_bin_s, assign last bin" << std::endl;
      std::cout << " last bin value: " << hVariable_Signal[i]     -> GetBinContent( bin_s ) << std::endl;
    }
    */

    if (bin_s>last_bin_s) bin_s = last_bin_s;
    if (bin_b>last_bin_b) bin_b = last_bin_b;

    bool kk= false;
    if (bin_s > last_bin_s){
      kk=true;
      std::cout << " the likelihood pdf value is  " << hVariable_Signal[i]->GetBinContent( bin_s ) << std::endl;
    }

    likelihood_s *= hVariable_Signal[i]     -> GetBinContent( bin_s );
    likelihood_b *= hVariable_Background[i] -> GetBinContent( bin_b );

    if (kk) std::cout << " likelihood s/b = " << likelihood_s << "/" << likelihood_b << std::endl;

  }

  float likelihood = (likelihood_s>0.) ? likelihood_s / ( likelihood_s + likelihood_b ) : 0.;

  // if (bin_s > last_bin_s)  std::cout << "final likelihood value = " << likelihood << std::endl;

  return likelihood;

}
