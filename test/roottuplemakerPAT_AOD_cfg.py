import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
############## IMPORTANT ########################################
# if you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100
#################################################################
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring('file:input_file.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RootTupleMakerPAT_output.root')
)

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Load additional RECO config
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V12::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# add tcMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')

## b-tagging general configuration
process.load("RecoBTag.Configuration.RecoBTag_cff")

## configure the softMuonByPt ESProducer and EDProducer
process.load("RecoBTag.SoftLepton.softLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.softMuonByPtBJetTags_cfi")

# b-tag discriminators to be added to the PAT jets
process.allLayer1Jets.discriminatorSources = cms.VInputTag(
    cms.InputTag("jetBProbabilityBJetTags"),
    cms.InputTag("simpleSecondaryVertexBJetTags"),
    cms.InputTag("softMuonByPtBJetTags"),
    cms.InputTag("trackCountingHighEffBJetTags"),
)

# Skimming
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

# RootTupleMaker
process.treeCreator = cms.EDAnalyzer('RootTupleMakerPAT')
process.treeCreator.maxgenparticles = cms.untracked.int32(25)
process.treeCreator.maxgenjets      = cms.untracked.int32(10)
process.treeCreator.maxelectrons    = cms.untracked.int32(10)
process.treeCreator.maxcalojets     = cms.untracked.int32(10)
process.treeCreator.maxmuons        = cms.untracked.int32(10)
process.treeCreator.debug           = cms.untracked.bool(False)
# overall luminosity normalization  (in pb-1)   
process.treeCreator.luminosity      = cms.untracked.double(100)
process.treeCreator.numEvents       = cms.untracked.int32(10)
process.treeCreator.saveTrigger     = cms.untracked.bool(True)
process.treeCreator.usePDFweight    = cms.untracked.bool(False)
process.treeCreator.PDFSet          = cms.untracked.string("/cteq61.LHgrid")
process.treeCreator.muonLabel       = cms.untracked.InputTag("cleanLayer1Muons");
process.treeCreator.electronLabel   = cms.untracked.InputTag("cleanLayer1Electrons");
process.treeCreator.caloJetLabel    = cms.untracked.InputTag("cleanLayer1Jets");
process.treeCreator.genJetLabel     = cms.untracked.InputTag("iterativeCone5GenJets");
process.treeCreator.electronIso     = cms.untracked.double(0.1);
process.treeCreator.muonIso         = cms.untracked.double(0.05);

# PAT sequence modification
process.patDefaultSequence.remove( process.countLayer1Objects )

# Switch off old trigger matching
from PhysicsTools.PatAlgos.tools.trigTools import switchOffTriggerMatchingOld
switchOffTriggerMatchingOld( process )

#process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path( process.softMuonByPtBJetTags*process.patDefaultSequence*process.LJFilter*process.treeCreator )

# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATTuple.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaDataForDroppedData = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('drop *')
)
process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands += patEventContent
process.out.outputCommands += [
    # PAT
    'keep *_layer1METs*_*_*',
    # GEN
    'keep recoGenParticles_genParticles_*_*',
    #'keep *_genEventScale_*_*',
    #'keep *_genEventWeight_*_*',
    'keep *_genEventPdfInfo_*_*',
    'keep *_genMet_*_*',
    'keep *_iterativeCone5GenJets_*_*',
    # TRIGGER
    'keep edmTriggerResults_TriggerResults_*_HLT',
    # PAT (dropped)
    'drop *_cleanLayer1Photons_*_*',
    'drop *_cleanLayer1Taus_*_*', 
    'drop *_cleanLayer1Hemispheres_*_*' ]

