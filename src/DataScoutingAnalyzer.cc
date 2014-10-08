// system include files
#include <memory>
#include <typeinfo>

// user include files
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "isildak/DataScoutingAnalyzer/interface/DataScoutingAnalyzer.h"

//objects
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "TLorentzVector.h"

template <typename jettype, typename mettype>
DataScoutingAnalyzer<jettype,mettype>::DataScoutingAnalyzer(const edm::ParameterSet& iConfig):
  tag_recoJet(iConfig.getParameter<edm::InputTag>("jets")),
  apply_corrections_reco(iConfig.getParameter<bool>("apply_corrections_reco")),
  apply_corrections_DS(iConfig.getParameter<bool>("apply_corrections_DS")),
  s_recoJetCorrector(iConfig.getParameter<std::string>("jetCorrectionsReco")),
  s_dsJetCorrector(iConfig.getParameter<std::string>("jetCorrectionsDS")),
  tag_recoRho(iConfig.getParameter<edm::InputTag>("rho")),  
  jetThreshold(iConfig.getParameter<double>("jetThreshold")),
  tag_recoMet(iConfig.getParameter<edm::InputTag>("met")),
  tag_recoElectrons(iConfig.getParameter<edm::InputTag>("electrons")),
  tag_recoMuons(iConfig.getParameter<edm::InputTag>("muons")),
  tag_hcalNoise(iConfig.getParameter<edm::InputTag>("noise")),
  s_outputFile(iConfig.getParameter<std::string>("outputFile")){

}

template <typename jettype, typename mettype>
DataScoutingAnalyzer<jettype,mettype>::~DataScoutingAnalyzer(){
}

// ------------ method called for each event  ------------
template <typename jettype, typename mettype>
void
DataScoutingAnalyzer<jettype,mettype>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  using namespace edm;

  runNo      = iEvent.id().run();
  evtNo      = iEvent.id().event();
  lumiBlock  = iEvent.id().luminosityBlock();

  Handle<std::vector<jettype> > h_recoJet;
  Handle<std::vector<mettype> > h_recoMet;
  const JetCorrector* correctorL1L2L3 = JetCorrector::getJetCorrector (s_recoJetCorrector, iSetup);

  const JetCorrector* correctorL2L3     = JetCorrector::getJetCorrector (s_dsJetCorrector, iSetup); 
 
  Handle<double> h_recoRho;
  edm::Handle<reco::JetIDValueMap> recoJetIDMap;

  iEvent.getByLabel(tag_recoJet,h_recoJet);
  iEvent.getByLabel(tag_recoMet,h_recoMet);
  iEvent.getByLabel(tag_recoRho,h_recoRho);

  Handle<std::vector<reco::CaloJet> > h_dsJet;
  Handle<std::vector<reco::CaloMET> > h_dsMet;
  Handle<std::vector<reco::CaloMET> > h_dsMetClean;
  Handle<double> h_dsRho;
 
  iEvent.getByLabel("hltCaloJetIDPassed",h_dsJet);
  iEvent.getByLabel("hltMet",h_dsMet);
  iEvent.getByLabel("hltMetClean",h_dsMetClean);
  iEvent.getByLabel("hltKT6CaloJets","rho",h_dsRho);
  
  edm::Handle< bool > HBHENoiseFilterResult;
  try { iEvent.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), HBHENoiseFilterResult); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << HBHENoiseFilterResult; }
  HBHENoiseFilterResultFlag = *HBHENoiseFilterResult;

  edm::Handle< bool > hcalLaserEventFilter;
  try { iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventFilter); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << hcalLaserEventFilter; }
  hcalLaserEventFilterFlag = *hcalLaserEventFilter;
  
  //fill the tree
  
  //MET
  dsMetPt = h_dsMet->front().pt();
  dsMetPhi = h_dsMet->front().phi();
  dsMetCleanPt = h_dsMetClean->front().pt();
  dsMetCleanPhi = h_dsMetClean->front().phi();
  recoMetPt = h_recoMet->front().pt();
  recoMetPhi = h_recoMet->front().phi();
  
  dsRho = *h_dsRho;
  recoRho = *h_recoRho;

  typename std::vector<jettype>::const_iterator i_recoJet;
  //std::vector<reco::Jet>::const_iterator i_recoJet;
  nRECOJets=0;
  for(i_recoJet = h_recoJet->begin(); i_recoJet != h_recoJet->end(); i_recoJet++)
  {  
    jettype recoJet = *i_recoJet;

    reco::CaloJet* p1 = dynamic_cast<reco::CaloJet*>(&recoJet);
    //reco::PFJet*   p2 = dynamic_cast<reco::PFJet*>(&recoJet);

    recoJetRawE[nRECOJets] = recoJet.energy();
    float tmp_pT = recoJet.pt();
    double scale = correctorL1L2L3->correction(recoJet,iEvent,iSetup);
    if(apply_corrections_reco)
    {
      recoJet.scaleEnergy(scale);
      recoJEC[nRECOJets] = recoJet.pt()/tmp_pT;
    }

    if(recoJet.pt() < jetThreshold) continue;
    
    recoJetPt[nRECOJets]    = recoJet.pt();
    recoJetEta[nRECOJets]   = recoJet.eta();
    recoJetPhi[nRECOJets]   = recoJet.phi();
    recoJetE[nRECOJets]     = recoJet.energy();

    if(p1)
    {
      recoJetFracHad[nRECOJets] = p1->energyFractionHadronic();
      recoJetFracEm[nRECOJets]  = p1->emEnergyFraction();
    }
    nRECOJets++;
  }

  std::vector<reco::CaloJet>::const_iterator i_dsJet;
  nDSJets=0;
  for(i_dsJet = h_dsJet->begin(); i_dsJet != h_dsJet->end(); i_dsJet++){
    reco::CaloJet dsJet = *i_dsJet;
    dsJetRawE[nDSJets]  = dsJet.energy();
    //apply pileup correction
    float pileupCorr      = 1-((1.08 - dsRho)*dsJet.jetArea()/dsJet.pt());

    if(dsJet.pt()*pileupCorr < jetThreshold) continue;

    if(apply_corrections_DS){

      if(pileupCorr > 0. && pileupCorr < 1.){
        dspileupCorr[nDSJets] = pileupCorr;
        dsJet.scaleEnergy(pileupCorr);
        float ds_tmp_pT = dsJet.pt();
        dsJet.scaleEnergy(correctorL2L3->correction(dsJet,iEvent,iSetup));
        dsJECL2L3Res[nDSJets] = dsJet.pt()/ds_tmp_pT;
      }
    }
    dsJetPt[nDSJets]      = dsJet.pt();
    dsJetEta[nDSJets]     = dsJet.eta();
    dsJetPhi[nDSJets]     = dsJet.phi();
    dsJetE[nDSJets]       = dsJet.energy();
    dsJetFracHad[nDSJets] = dsJet.energyFractionHadronic();
    dsJetFracEm[nDSJets]  = dsJet.emEnergyFraction();

    //do jet matching
    float bestdEoE  = 9999;
    int bestIndex   = -1;
    float minDeltaR = 9999.;
    float dR;
    for( int iRECOJet=0; iRECOJet < nRECOJets; iRECOJet++){
      dR= reco::deltaR(dsJet.eta(),dsJet.phi(),recoJetEta[iRECOJet],recoJetPhi[iRECOJet]);
      if( dR> 0.5) continue; //require DR match
      if (dR < minDeltaR){
        minDeltaR = dR;
        bestIndex = iRECOJet;
      }

      float dEoE = fabs(dsJet.energy() - recoJetE[iRECOJet])/recoJetE[iRECOJet];

      if(dEoE < 0.5 && dEoE < bestdEoE){
        bestdEoE = dEoE;
        bestIndex = iRECOJet;
      }
    }
    dsJetMatchIndex[nDSJets] = bestIndex;
    nDSJets++;
  }
  outputTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::beginJob(){

  outputFile = new TFile(s_outputFile.c_str(),"RECREATE");
  outputTree = new TTree("DSComp","");

  outputTree->Branch("runNo",&runNo,"runNo/I");
  outputTree->Branch("evtNo",&evtNo,"evtNo/I");
  outputTree->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");

  outputTree->Branch("nDSJets",&nDSJets,"nDSJets/I");

  outputTree->Branch("dspileupCorr",dspileupCorr,"dspileupCorr[nDSJets]");
  outputTree->Branch("recoJEC",recoJEC,"recoJEC[nDSJets]");
  outputTree->Branch("dsJECL2L3Res",dsJECL2L3Res,"dsJECL2L3Res[nDSJets]");
  outputTree->Branch("dsJetPt",dsJetPt,"dsJetPt[nDSJets]");
  outputTree->Branch("dsJetEta",dsJetEta,"dsJetEta[nDSJets]");
  outputTree->Branch("dsJetPhi",dsJetPhi,"dsJetPhi[nDSJets]");
  outputTree->Branch("dsJetE",dsJetE,"dsJetE[nDSJets]");
  outputTree->Branch("dsJetRawE",dsJetRawE,"dsJetRawE[nDSJets]");
  outputTree->Branch("dsJetFracHad",dsJetFracHad,"dsJetFracHad[nDSJets]");
  outputTree->Branch("dsJetFracEm",dsJetFracEm,"dsJetFracEm[nDSJets]");
  outputTree->Branch("dsJetMatchIndex",dsJetMatchIndex,"dsJetMatchIndex[nDSJets]/I");

  outputTree->Branch("dsRho",&dsRho);
  outputTree->Branch("dsMetPt",&dsMetPt);
  outputTree->Branch("dsMetPhi",&dsMetPhi);
  outputTree->Branch("dsMetCleanPt",&dsMetCleanPt);
  outputTree->Branch("dsMetCleanPhi",&dsMetCleanPhi);

  outputTree->Branch("nRECOJets",&nRECOJets,"nRECOJets/I");
  outputTree->Branch("recoJetPt",recoJetPt,"recoJetPt[nRECOJets]");
  outputTree->Branch("recoJetEta",recoJetEta,"recoJetEta[nRECOJets]");
  outputTree->Branch("recoJetPhi",recoJetPhi,"recoJetPhi[nRECOJets]");
  outputTree->Branch("recoJetE",recoJetE,"recoJetE[nRECOJets]");
  outputTree->Branch("recoJetRawE",recoJetRawE,"recoJetRawE[nRECOJets]");
  outputTree->Branch("recoJetE",recoJetE,"recoJetE[nRECOJets]");
  outputTree->Branch("recoJetFracHad",recoJetFracHad,"recoJetFracHad[nRECOJets]");
  outputTree->Branch("recoJetFracEm",recoJetFracEm,"recoJetFracEm[nRECOJets]");

  outputTree->Branch("recoRho",&recoRho);
  outputTree->Branch("recoMetPt",&recoMetPt);
  outputTree->Branch("recoMetPhi",&recoMetPhi);
  outputTree->Branch("recoMetCleanPt",&recoMetCleanPt);
  outputTree->Branch("recoMetCleanPhi",&recoMetCleanPhi);

  outputTree->Branch("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlag,"HBHENoiseFilterResultFlag/B");
  outputTree->Branch("eeBadScFilterFlag",&eeBadScFilterFlag,"eeBadScFilterFlag/B"); 
  outputTree->Branch("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlag,"hcalLaserEventFilterFlag/B"); 
}

// ------------ method called once each job just after ending the event loop  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::endJob() 
{
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
}

// ------------ method called when starting to processes a run  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::beginRun(edm::Run const&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a run  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::endRun(edm::Run const&, edm::EventSetup const&){
}

// ------------ method called when starting to processes a luminosity block  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a luminosity block  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <typename jettype, typename mettype>
void
DataScoutingAnalyzer<jettype,mettype>::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


typedef DataScoutingAnalyzer<reco::CaloJet,reco::CaloMET> CaloScoutingAnalyzer;
typedef DataScoutingAnalyzer<reco::PFJet,reco::CaloMET> PFJetScoutingAnalyzer;
typedef DataScoutingAnalyzer<reco::PFJet,reco::PFMET> PFScoutingAnalyzer;
//define this as a plug-in
DEFINE_FWK_MODULE(CaloScoutingAnalyzer);
DEFINE_FWK_MODULE(PFJetScoutingAnalyzer);
DEFINE_FWK_MODULE(PFScoutingAnalyzer);
