import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")

process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.taggingMode = True
process.hcalLaserEventFilter.vetoByRunEventNumber = False
process.hcalLaserEventFilter.vetoByHBHEOccupancy= True

process.load("RecoMET.METFilters.eeBadScFilter_cfi")
process.eeBadScFilter.taggingMode = True



process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_0.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_1.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_10.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_11.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_12.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_13.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_14.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_15.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_16.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_17.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_18.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_19.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_2.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_20.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_21.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_22.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_23.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_24.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_25.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_26.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_27.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_28.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_29.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_3.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_30.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_31.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_32.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_33.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_34.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_35.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_36.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_37.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_38.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_39.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_4.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_40.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_41.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_42.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_43.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_44.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_45.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_46.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_47.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_5.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_6.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_7.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_8.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_9.root',
)
)

VERTEX_SEL=("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")

process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
                                   src = cms.InputTag("offlinePrimaryVertices"),
                                   cut = cms.string(VERTEX_SEL),
                                   filter = cms.bool(True),
                                   )

plugin_type       = [["CaloScoutingAnalyzer","ak5CaloJets","ak5CaloL1L2L3Residual","kt6CaloJets"],["PFJetScoutingAnalyzer","ak5PFJets","ak5PFL1FastL2L3Residual","kt6PFJets"]]
i = 0 # 0 is CaloScoutingAnalyzer 1 is PFJetScoutingAnalyzer

print plugin_type[i][0]+" is running."+"\njets:"+plugin_type[i][1]+"\njetCorrections:"+plugin_type[i][2]+"\nrho:"+plugin_type[i][3]
if i==0:
	output = 'test_Calo.root'
else:
	output = 'test_PF.root'
process.demo = cms.EDAnalyzer(plugin_type[i][0],
                              jets = cms.InputTag(plugin_type[i][1]),
                              HBHENoiseFilterResultLabel = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"),
                              apply_corrections_reco = cms.bool(True),
                              apply_corrections_DS = cms.bool(True),
                              jetCorrectionsReco = cms.string(plugin_type[i][2]),
                              jetCorrectionsDS = cms.string("ak5CaloL2L3Residual"),
                              rho  = cms.InputTag(plugin_type[i][3],"rho"),
                              jetThreshold = cms.double(20),
                              met = cms.InputTag("met"),
                              electrons = cms.InputTag(""),
                              muons = cms.InputTag(""),
                              noise = cms.InputTag(""),
                              outputFile = cms.string(output)

)

process.metOptionalFilterSequence = cms.Sequence(process.HBHENoiseFilterResultProducer * process.hcalLaserEventFilter * process.eeBadScFilter)
#process.goodPrimaryVertices *
process.p       = cms.Path(  process.goodPrimaryVertices * process.metOptionalFilterSequence * process.demo)

