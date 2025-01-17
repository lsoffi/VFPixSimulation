import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.TrackJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process = cms.Process("TrkJetPhase2")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/TP2023SHCALDR/VBF_HToInv_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/40000/0485645C-62DB-E411-A8C1-0025905B8562.root'
    )
)


process.ak4TrackJets = cms.EDProducer(
    "FastjetJetProducer",
    TrackJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4)
    )




process.TrkJetPhase2 = cms.EDAnalyzer('TrkJetAnalyzer',
                                      GenJetCollection       = cms.InputTag("ak5GenJets"),
                                      GenParticlesCollection = cms.InputTag("genParticles"),
                                      TrackJetCollection    = cms.InputTag("ak4TrackJets"),
                                      TrackCollection    = cms.InputTag("generalTracks"),
                                     VertexCollection = cms.InputTag("offlinePrimaryVertices"),
                                      PUSummaryInfoCollection = cms.InputTag("addPileupInfo"),
                                      minTracks = cms.untracked.uint32(1000)
         )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histoOUTPUT.root')
                                   )
#dump event content
process.dump=cms.EDAnalyzer('EventContentAnalyzer')


process.p = cms.Path(process.ak4TrackJets*process.TrkJetPhase2)
