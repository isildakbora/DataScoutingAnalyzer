#ifndef DSComp_h
#define DSComp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class DSComp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runNo;
   Int_t           evtNo;
   Int_t           lumiBlock;
   Int_t           nDSJets;
   Float_t         dsJetPt[35];   //[nDSJets]
   Float_t         dsJetEta[35];   //[nDSJets]
   Float_t         dsJetPhi[35];   //[nDSJets]
   Float_t         dsJetE[35];   //[nDSJets]
   Float_t         dsJetRawE[35];   //[nDSJets]
   Float_t         dspileupCorr[35];   //[nDSJets]
   Float_t         dsJECL2L3Res[35];   //[nDSJets]
   Float_t         dsJetFracHad[35];   //[nDSJets]
   Float_t         dsJetFracEm[35];   //[nDSJets]
   Int_t           dsJetMatchIndex[35];   //[nDSJets]
   Float_t         dsRho;
   Float_t         dsMetPt;
   Float_t         dsMetPhi;
   Float_t         dsMetCleanPt;
   Float_t         dsMetCleanPhi;
   Int_t           nRECOJets;
   Float_t         recoJetPt[89];   //[nRECOJets]
   Float_t         recoJEC[89];   //[nRECOJets]
   Float_t         recoJetEta[89];   //[nRECOJets]
   Float_t         recoJetPhi[89];   //[nRECOJets]
   Float_t         recoJetE[89];   //[nRECOJets]
   Float_t         recoJetRawE[89];   //[nRECOJets]
   Float_t         recoRho;
   Float_t         recoMetPt;
   Float_t         recoMetPhi;
   bool            HBHENoiseFilterResultFlag;
   bool            eeBadScFilterFlag;
   bool            hcalLaserEventFilterFlag; 

   // List of branches
   TBranch        *b_runNo;   //!
   TBranch        *b_evtNo;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_nDSJets;   //!
   TBranch        *b_dspileupCorr; //!
   TBranch        *b_recoJEC;//!
   TBranch        *b_dsJECL2L3Res;//!
   TBranch        *b_dsJetPt;   //!
   TBranch        *b_dsJetEta;   //!
   TBranch        *b_dsJetPhi;   //!
   TBranch        *b_dsJetE;   //!
   TBranch        *b_dsJetRawE;   //!
   TBranch        *b_dsJetFracHad;   //!
   TBranch        *b_dsJetFracEm;   //!
   TBranch        *b_dsJetMatchIndex;   //!
   TBranch        *b_dsRho;   //!
   TBranch        *b_dsMetPt;   //!
   TBranch        *b_dsMetPhi;   //!
   TBranch        *b_dsMetCleanPt;   //!
   TBranch        *b_dsMetCleanPhi;   //!
   TBranch        *b_nRECOJets;   //!
   TBranch        *b_recoJetPt;   //!
   TBranch        *b_recoJetEta;   //!
   TBranch        *b_recoJetPhi;   //!
   TBranch        *b_recoJetE;   //!
   TBranch        *b_recoJetRawE;   //!
   TBranch        *b_recoJetE;   //!
   TBranch        *b_recoRho;   //!
   TBranch        *b_recoMetPt;   //!
   TBranch        *b_recoMetPhi;   //!
   TBranch        *b_HBHENoiseFilterResultFlag; //!
   TBranch        *b_eeBadScFilterFlag; //!
   TBranch        *b_hcalLaserEventFilterFlag; //!

   DSComp(TTree *tree=0);
   virtual ~DSComp();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DSComp_cxx
DSComp::DSComp(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TString filename[2]= {"test_Calo.root","test_PF.root"};
   int typeno=0;
   cout << filename[typeno] << endl;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename[typeno]);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename[typeno]);
      }
      f->GetObject("DSComp",tree);

   }
   Init(tree);
}

DSComp::~DSComp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DSComp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DSComp::LoadTree(Long64_t entry)
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

void DSComp::Init(TTree *tree)
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

   fChain->SetBranchAddress("runNo", &runNo, &b_runNo);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evtNo);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("nDSJets", &nDSJets, &b_nDSJets);
   fChain->SetBranchAddress("dspileupCorr", dspileupCorr, &b_dspileupCorr);
   fChain->SetBranchAddress("recoJEC", recoJEC, b_recoJEC);
   fChain->SetBranchAddress("dsJECL2L3Res", dsJECL2L3Res, b_dsJECL2L3Res);
   fChain->SetBranchAddress("dsJetPt", dsJetPt, &b_dsJetPt);
   fChain->SetBranchAddress("dsJetEta", dsJetEta, &b_dsJetEta);
   fChain->SetBranchAddress("dsJetPhi", dsJetPhi, &b_dsJetPhi);
   fChain->SetBranchAddress("dsJetE", dsJetE, &b_dsJetE);
   fChain->SetBranchAddress("dsJetRawE", dsJetRawE, &b_dsJetRawE);
   fChain->SetBranchAddress("dsJetFracHad", dsJetFracHad, &b_dsJetFracHad);
   fChain->SetBranchAddress("dsJetFracEm", dsJetFracEm, &b_dsJetFracEm);
   fChain->SetBranchAddress("dsJetMatchIndex", dsJetMatchIndex, &b_dsJetMatchIndex);
   fChain->SetBranchAddress("dsRho", &dsRho, &b_dsRho);
   fChain->SetBranchAddress("dsMetPt", &dsMetPt, &b_dsMetPt);
   fChain->SetBranchAddress("dsMetPhi", &dsMetPhi, &b_dsMetPhi);
   fChain->SetBranchAddress("dsMetCleanPt", &dsMetCleanPt, &b_dsMetCleanPt);
   fChain->SetBranchAddress("dsMetCleanPhi", &dsMetCleanPhi, &b_dsMetCleanPhi);
   fChain->SetBranchAddress("nRECOJets", &nRECOJets, &b_nRECOJets);
   fChain->SetBranchAddress("recoJetPt", recoJetPt, &b_recoJetPt);
   fChain->SetBranchAddress("recoJetEta", recoJetEta, &b_recoJetEta);
   fChain->SetBranchAddress("recoJetPhi", recoJetPhi, &b_recoJetPhi);
   fChain->SetBranchAddress("recoJetE", recoJetE, &b_recoJetE);
   fChain->SetBranchAddress("recoJetRawE", recoJetRawE, &b_recoJetRawE);
   fChain->SetBranchAddress("recoRho", &recoRho, &b_recoRho);
   fChain->SetBranchAddress("recoMetPt", &recoMetPt, &b_recoMetPt);
   fChain->SetBranchAddress("recoMetPhi", &recoMetPhi, &b_recoMetPhi);
   fChain->SetBranchAddress("HBHENoiseFilterResultFlag", &HBHENoiseFilterResultFlag, &b_HBHENoiseFilterResultFlag);
   fChain->SetBranchAddress("eeBadScFilterFlag", &eeBadScFilterFlag, &b_eeBadScFilterFlag);
   fChain->SetBranchAddress("hcalLaserEventFilterFlag", &hcalLaserEventFilterFlag, &b_hcalLaserEventFilterFlag);
   Notify();
}

Bool_t DSComp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DSComp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DSComp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DSComp_cxx
