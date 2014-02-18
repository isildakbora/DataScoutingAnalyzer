#define DSComp_cxx
#include "DSComp.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <string>
void DSComp::Loop()
{
  //create a tree file 
  TFile *mjjtree = new TFile("dijet.root","recreate");

  TTree *DijetMassTree = new TTree("DijetMassTree","DijetMassTree");

  TFile *f = new TFile("Delta_pT_histograms.root", "RECREATE");

  TH1F* h_recoJetPt      = new TH1F("RecoJetpT","Reco Jet pT",498,20,5000);

  TH1F* h_recoJetE       = new TH1F("RecoJetE","Reco Jet E",498,20,5000);
  TH1F* h_dsJetE         = new TH1F("HLTJetE","HLT Jet E",498,20,5000);

  TH1F* h_recoJetRawE    = new TH1F("RecoJetRawE","Reco Jet Raw E",498,20,5000);
  TH1F* h_dsJetRawE      = new TH1F("HLTJetRawE","HLT Jet Raw E",498,20,5000);

  TH1F* h_recoHLTpTdiff  = new TH1F("recoHLTpTdiff","recoHLTpTdiff",400,-1.,1.);
  TH1F* h_recoHLTEtadiff = new TH1F("recoHLTEtadiff","recoHLTEtadiff",4800,-2,2);
  TH1F* h_recoHLTPhidiff = new TH1F("recoHLTPhidiff","recoHLTPhidiff",4800,-TMath::Pi(), TMath::Pi());


  Float_t recoMjj, dsMjj, dsJetPt0, dsJetE0, dsJetE1, dsJetPt1, dsJetEta0, dsJetEta1, dsJetPhi0, dsJetPhi1,
   recoJetPt0, recoJetPt1, recoJetE0, recoJetE1, recoJetEta0, recoJetEta1, recoJetPhi0, recoJetPhi1;

  DijetMassTree->Branch("runNo",&runNo,"runNo/I");   
  DijetMassTree->Branch("evtNo",&evtNo,"evtNo/I");
  DijetMassTree->Branch("dsMjj",&dsMjj,"dsMjj/F");
  DijetMassTree->Branch("recoMjj",&recoMjj,"recoMjj/F");


  DijetMassTree->Branch("dsJetEta0",&dsJetEta0,"dsJetEta0/F");
  DijetMassTree->Branch("dsJetEta1",&dsJetEta1,"dsJetEta1/F");

  DijetMassTree->Branch("dsJetPhi0",&dsJetPhi0,"dsJetPhi0/F");
  DijetMassTree->Branch("dsJetPhi1",&dsJetPhi1,"dsJetPhi1/F");

  DijetMassTree->Branch("dsJetPt0",&dsJetPt0,"dsJetPt0/F");
  DijetMassTree->Branch("dsJetPt1",&dsJetPt1,"dsJetPt1/F");

  DijetMassTree->Branch("dsJetE0",&dsJetE0,"dsJetE0/F");
  DijetMassTree->Branch("dsJetE1",&dsJetE1,"dsJetE1/F");

  DijetMassTree->Branch("recoJetEta0",&recoJetEta0,"recoJetEta0/F");
  DijetMassTree->Branch("recoJetEta1",&recoJetEta1,"recoJetEta1/F");
  
  DijetMassTree->Branch("recoJetPhi0",&recoJetEta0,"recoJetPhi0/F");
  DijetMassTree->Branch("recoJetPhi1",&recoJetPhi1,"recoJetPhi1/F");

  DijetMassTree->Branch("recoJetPt0",&recoJetPt0,"recoJetPt0/F");
  DijetMassTree->Branch("recoJetPt1",&recoJetPt1,"recoJetPt1/F");

  DijetMassTree->Branch("recoJetE0",&recoJetE0,"recoJetE0/F");
  DijetMassTree->Branch("recoJetE1",&recoJetE1,"recoJetE1/F");

  char *HLT[2]  = {"HLT Pass","HLT Fail"};
  char *RECO[2] = {"RECO Pass","RECO Fail"};

  TH2F* h_HLT_RECO_Table = new TH2F("HLT_RECO_Table", "HLT_RECO_Table",2,0,2,2,0,2);
  for (int i=1;i<=2;i++) h_HLT_RECO_Table->GetXaxis()->SetBinLabel(i,HLT[i-1]);
  for (int i=1;i<=2;i++) h_HLT_RECO_Table->GetYaxis()->SetBinLabel(i,RECO[i-1]);

   if (fChain == 0) return;

   Float_t dsjetpt, recojetpt, dsjeteta, recojeteta;
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << nentries << std::endl;
   Long64_t nbytes = 0, nb = 0;
    
    

   // Set kinematic binnings

    float pTbins[]  = {30., 103., 206., 309., 413., 515., 618., 721., 1237.};
    const Int_t   n_pTbins  = sizeof(pTbins)/sizeof(float)-1;

    float etaBins[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5};
    const Int_t   n_etaBins = sizeof(etaBins)/sizeof(float)-1;
    TString histoname;
    TH1F *hDeltapT[8][5];
    TH1F *hDeltaMjj[8][5];
    TH1F *hDeltapT12[8][5];
    TH1F *hdspileupCorr[8][5];
    TH1F *hdsJECL2L3Res[8][5];
    TH1F *hrecoJEC[8][5];
    TF1 *cbfit[8][5];
    TF1 *aux_xBal[8][5];
    TF1 *aux_gaus[8][5];

  // Set the names and bins of histograms
    for (int i = 0; i < n_pTbins; ++i)
    {
        for (int j=0; j<n_etaBins; ++j)
        {
          histoname = "Delta_pT_"+TString::Format("%.0f",pTbins[i])+"_"+TString::Format("%.0f",pTbins[i+1])+"_eta_"+TString::Format("%2.1f",j*0.5)+"_"+TString::Format("%2.1f",(j+1)*0.5);
          hDeltapT[i][j]= new TH1F(histoname, histoname, 1200, -1.2, 1.2);
          hDeltapT[i][j]->Sumw2();

          hDeltaMjj[i][j]= new TH1F(histoname+"Mjj", histoname+"Mjj", 1200, -1.2, 1.2);
          hDeltaMjj[i][j]->Sumw2();

          hDeltapT12[i][j]= new TH1F(histoname+"pT12", histoname+"pT12", 1200, -1.2, 1.2);
          hDeltapT12[i][j]->Sumw2();

          histoname = "dspileopCorr"+TString::Format("%.0f",pTbins[i])+"_"+TString::Format("%.0f",pTbins[i+1])+"_eta_"+TString::Format("%2.1f",j*0.5)+"_"+TString::Format("%2.1f",(j+1)*0.5);
          hdspileupCorr[i][j]= new TH1F(histoname, histoname, 1000, 0, 10);

          histoname = "dsJECL2L3Res"+TString::Format("%.0f",pTbins[i])+"_"+TString::Format("%.0f",pTbins[i+1])+"_eta_"+TString::Format("%2.1f",j*0.5)+"_"+TString::Format("%2.1f",(j+1)*0.5);
          hdsJECL2L3Res[i][j]= new TH1F(histoname, histoname, 1000, 0, 10);

          histoname = "recoJEC"+TString::Format("%.0f",pTbins[i])+"_"+TString::Format("%.0f",pTbins[i+1])+"_eta_"+TString::Format("%2.1f",j*0.5)+"_"+TString::Format("%2.1f",(j+1)*0.5);
          hrecoJEC[i][j]= new TH1F(histoname, histoname, 1000, 0, 10);
        }  
    }

   // Set threshold values at which events are excluded
    float minPtThreshold = 30.; // GeV
    float maxAbsEtaThreshold = 2.5;
    float maxEtaSepThreshold = 2.0;
    
   // Tree Loop
   int decade = 0;
   double progress = 0;
   for (Long64_t jentry=0; jentry<nentries/1; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //----------- progress report -------------
      double progress = 100.*jentry/nentries;
      int k = TMath::FloorNint(progress);
      std::cout << "\r" << progress << "% completed: ";
      //std::cout << "[" << std::string(k, '|') << std::string(99-k, ' ') << "]";
      std::cout.flush();
      decade = k; 

      // Filters //
      //if two jets exist
      if (!((nDSJets > 1) && (nRECOJets > 1))) continue;    

      //Dijet HLT Selection
      int DeltaEtaHLT           = (fabs(dsJetEta[0]-dsJetEta[1]) < 1.3) ? 1 : 0;
      int MaxAbsEtaThresholdHLT = ((fabs(dsJetEta[0]) < maxAbsEtaThreshold) && (fabs(dsJetEta[1]) < maxAbsEtaThreshold));
      int MinPtThresholdHLT     = ((dsJetPt[0] > minPtThreshold) && (dsJetPt[1] > minPtThreshold));
      int allHLTDijetSelection  = (DeltaEtaHLT && MinPtThresholdHLT && MaxAbsEtaThresholdHLT);     

      //Dijet RECO Selection
      int matchindex0 = dsJetMatchIndex[0];
      int matchindex1 = dsJetMatchIndex[1];
      int recoJetID;
      int allRECODijetSelection;
      if(matchindex0 >=0 && matchindex1 >=0)
      {
        int DeltaEtaRECO           = (fabs(recoJetEta[matchindex0]-recoJetEta[matchindex1]) < 1.3) ? 1 : 0;
        int MinPtThresholdRECO     = ((recoJetPt[matchindex0] > minPtThreshold) && (recoJetPt[matchindex1] > minPtThreshold));
        int MaxAbsEtaThresholdRECO = ((fabs(recoJetEta[matchindex0]) < maxAbsEtaThreshold) && (fabs(recoJetEta[matchindex1]) < maxAbsEtaThreshold));
        int recoJetID              = (recoJetFracHad[matchindex0] < 0.95 && recoJetFracHad[matchindex1] < 0.95 && recoJetFracEm[matchindex0] < 0.95 && recoJetFracEm[matchindex1] < 0.95) ? 1 : 0;
        allRECODijetSelection      = (DeltaEtaRECO && MinPtThresholdRECO && MaxAbsEtaThresholdRECO);  
      }
      else  allRECODijetSelection = 0;
      //RECO Event Filters
      int DeltaPhiRECO  = (fabs(recoJetPhi[0]-recoJetPhi[1]) > TMath::Pi()/3) ? 1 : 0;
      int RecoFlagsGood = (DeltaPhiRECO && HBHENoiseFilterResultFlag && hcalLaserEventFilterFlag && eeBadScFilterFlag);

      //HLT Event Filters
      int DeltaPhiHLT         = (fabs(dsJetPhi[0]-dsJetPhi[1]) > TMath::Pi()/3) ? 1 : 0;
      int MET_vs_METCleanFlag = (dsMetPt != dsMetCleanPt) ? 0 : 1;
      int dsJetID             = (dsJetFracHad[0] < 0.95 && dsJetFracHad[1] < 0.95 && dsJetFracEm[0] < 0.95 && dsJetFracEm[1] < 0.95) ? 1 : 0;
      int HLTFlagsGood        = (dsJetID && MET_vs_METCleanFlag && DeltaPhiHLT);

      int keepEvent = (allHLTDijetSelection && allRECODijetSelection);

      if (!keepEvent)
      {
        continue; 
      }

      // Fill the HLT vs RECO Filter Table
      if(HLTFlagsGood && RecoFlagsGood)        h_HLT_RECO_Table->Fill(0,0);
      else if(!HLTFlagsGood && !RecoFlagsGood) h_HLT_RECO_Table->Fill(1,1);
      else if(!HLTFlagsGood && RecoFlagsGood)  h_HLT_RECO_Table->Fill(1,0);
      else                                     h_HLT_RECO_Table->Fill(0,1);

        // Filters //

      // Dijet mass difference
      int pTbin, etabin;
      double dsX0, dsY0, dsZ0, dsX1, dsY1, dsZ1,
      recoX0, recoY0, recoZ0, recoX1, recoY1, recoZ1;
      if(matchindex0 >=0 && matchindex1 >=0)
      {
        // ds Jets
        dsJetPt0  = dsJetPt[0];
        dsJetPt1  = dsJetPt[1];
        dsJetEta0 = dsJetEta[0];
        dsJetEta1 = dsJetEta[1];
        dsJetPhi0 = dsJetPhi[0];
        dsJetPhi1 = dsJetPhi[1];

        dsJetE0 = dsJetE[0];                           
        dsJetE1 = dsJetE[1];

        dsX0 = dsJetPt[0]*TMath::Cos(dsJetPhi0);
        dsX1 = dsJetPt[1]*TMath::Cos(dsJetPhi1);
        dsY0 = dsJetPt[0]*TMath::Sin(dsJetPhi0);
        dsY1 = dsJetPt[1]*TMath::Sin(dsJetPhi1);
        dsZ0 = dsJetPt[0]*TMath::SinH(dsJetEta0);
        dsZ1 = dsJetPt[1]*TMath::SinH(dsJetEta1);

        // RECO Jets
        recoJetPt0  = recoJetPt[matchindex0];
        recoJetPt1  = recoJetPt[matchindex1];
        recoJetEta0 = recoJetEta[matchindex0];
        recoJetEta1 = recoJetEta[matchindex1];
        recoJetPhi0 = recoJetPhi[matchindex0];
        recoJetPhi1 = recoJetPhi[matchindex1];

        recoJetE0 = recoJetE[matchindex0];
        recoJetE1 = recoJetE[matchindex1];

        recoX0 = recoJetPt0*TMath::Cos(recoJetPhi0);
        recoX1 = recoJetPt1*TMath::Cos(recoJetPhi1);
        recoY0 = recoJetPt0*TMath::Sin(recoJetPhi0);
        recoY1 = recoJetPt1*TMath::Sin(recoJetPhi1);
        recoZ0 = recoJetPt0*TMath::SinH(recoJetEta0);
        recoZ1 = recoJetPt1*TMath::SinH(recoJetEta1);

        dsMjj   = TMath::Sqrt((dsJetE0+dsJetE1)**2 - (dsX0+dsX1)**2 - (dsY0+dsY1)**2 - (dsZ0+dsZ1)**2);
        recoMjj = TMath::Sqrt((recoJetE0+recoJetE1)**2 - (recoX0+recoX1)**2 - (recoY0+recoY1)**2 - (recoZ0+recoZ1)**2);

        pTbin     = Get_ij(pTbins, n_pTbins, etaBins, n_etaBins, dsPt,dsEta).first;
        etabin    = Get_ij(pTbins, n_pTbins, etaBins, n_etaBins, dsPt,dsEta).second;
        hDeltaMjj[pTbin][etabin]->Fill((dsMjj-recoMjj)/recoMjj);
        DijetMassTree->Fill();

        hDeltapT12[pTbin][etabin]->Fill((dsJetPt0-recoJetPt0)/recoJetPt0);
        hDeltapT12[pTbin][etabin]->Fill((dsJetPt1-recoJetPt0)/recoJetPt1);
      }

  // Jet Loop//
  for(Int_t i=0; i< nDSJets; i++)
  {
     Double_t frac_diff, dsPt, recoPt, dsEta, recoEta, dsPhi, recoPhi;
     int matchindex;
     
      if( nDSJets<35)
      {
        matchindex = dsJetMatchIndex[i];
        if( matchindex >= 0)
        {
          dsPt      = dsJetPt[i];
          recoPt    = recoJetPt[matchindex];
          dsEta     = dsJetEta[i];
          recoEta   = recoJetEta[matchindex];
          dsPhi     = dsJetPhi[i];
          recoPhi   = recoJetPhi[matchindex];

          frac_diff = (dsPt-recoPt)/recoPt;

          h_dsJetRawE   -> Fill(dsJetRawE[i]);
          h_recoJetRawE -> Fill(recoJetRawE[i]);

          h_dsJetE      -> Fill(dsJetE[i]);
          h_recoJetE    -> Fill(recoJetE[i]);

          h_recoHLTEtadiff -> Fill((dsEta - recoEta)/recoEta);
          h_recoHLTPhidiff -> Fill((dsPhi - recoPhi)/recoPhi);
          h_recoHLTpTdiff  -> Fill(frac_diff);


          pTbin     = Get_ij(pTbins, n_pTbins, etaBins, n_etaBins, dsPt,dsEta).first;
          etabin    = Get_ij(pTbins, n_pTbins, etaBins, n_etaBins, dsPt,dsEta).second;

          hDeltapT[pTbin][etabin]->Fill(frac_diff);

          hdspileupCorr[pTbin][etabin] -> Fill(dspileupCorr[i]);
          hdsJECL2L3Res[pTbin][etabin] -> Fill(dsJECL2L3Res[i]);
          hrecoJEC[pTbin][etabin]      -> Fill(recoJEC[i]);
        }
      }
  }
        // Jet Loop//
  }
   // Tree Loop

  TString name;
  f->cd();
  Double_t maxVal, mean, RMS;
  for (int i = 0; i < n_pTbins; ++i)
  {
    for (int j=0; j<n_etaBins; ++j)
    {
      name = "Delta_pT_"+TString::Format("%.0f",pTbins[i])+"_"+TString::Format("%.0f",pTbins[i+1])+"_eta_"+TString::Format("%2.1f",j*0.5)+"_"+TString::Format("%2.1f",(j+1)*0.5);
            
      hDeltapT[i][j]  ->Scale(1./hDeltapT[i][j]->GetEntries());
      hDeltaMjj[i][j]->Scale(1./hDeltaMjj[i][j]->GetEntries());
            
      hDeltapT[i][j]      ->Write();
      hDeltaMjj[i][j]    ->Write();
      hdspileupCorr[i][j] ->Write();
      hdsJECL2L3Res[i][j] ->Write();
      hrecoJEC[i][j]      ->Write();
    }  
           
  }
  h_HLT_RECO_Table ->Write();
  h_dsJetRawE   -> Write();
  h_recoJetRawE -> Write();
  h_dsJetE      -> Write();
  h_recoJetE    -> Write();
  f-> Write();  
  mjjtree->cd();
  DijetMassTree->Write();
  mjjtree->Close();
}

std::pair <int,int> Get_ij(float array1[], int length1, float array2[], int length2, double num1, double num2)
{
    int i = 0;
    int j = 0;
    
    for(int k=0; k < length1; ++k)
    {
        if(num1 > array1[k])
        {
            i = k;
        }
        else if(num1 > array1[k] && num1 < array1[k+1])
            break;
        
    }
    
    for(int k=0; k < length2; ++k)
    {
        if(num2 > array2[k])
        {
            j = k;
        }
        else if(num2 > array2[k] && num1 < array2[k+1])
            break;
    }
    return std::make_pair(i,j);
}
