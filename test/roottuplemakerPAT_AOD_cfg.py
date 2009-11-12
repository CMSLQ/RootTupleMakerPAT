import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
############## IMPORTANT ########################################
# if you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 10
#################################################################
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Source
process.source = cms.Source("PoolSource",
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck'), 
     fileNames = cms.untracked.vstring('file:input_file.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RootTupleMakerPAT_output.root')
)

# Load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# HEEPify the PAT electrons
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
    eleLabel = cms.InputTag("allLayer1Electrons"),
    barrelCuts = cms.PSet(heepBarrelCuts),
    endcapCuts = cms.PSet(heepEndcapCuts)
)

# Add 'heepPatElectrons' in the right place and point 'selectedLayer1Electrons' to them
process.patDefaultSequence.replace( process.allLayer1Electrons, process.allLayer1Electrons*process.heepPatElectrons )
process.selectedLayer1Electrons.src = cms.InputTag("heepPatElectrons")

# Electron and jet cleaning deltaR parameters
process.cleanLayer1Electrons.checkOverlaps.muons.deltaR = 0.3
process.cleanLayer1Jets.checkOverlaps.muons.deltaR = 0.5
process.cleanLayer1Jets.checkOverlaps.electrons.deltaR = 0.5

# Load additional RECO config
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All'
process.load("Configuration.StandardSequences.MagneticField_cff")

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# b-tag discriminators to be added to the PAT jets
process.allLayer1Jets.discriminatorSources = cms.VInputTag(
    cms.InputTag("jetBProbabilityBJetTags"),
    cms.InputTag("simpleSecondaryVertexBJetTags"),
    cms.InputTag("softMuonByPtBJetTags"),
    cms.InputTag("trackCountingHighEffBJetTags"),
)

# Fix InputTags for ECAL IsoDeposits (to work with FastSim samples)
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")
process.gamIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
process.gamIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")

# Load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.coreTools import *
restrictInputToAOD(process, ['All'])

# Skim definition
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
process.LJFilter.muLabel = 'muons'
process.LJFilter.elecLabel = 'gsfElectrons'
process.LJFilter.jetLabel = 'iterativeCone5CaloJets'
process.LJFilter.muPT = 20.
process.LJFilter.elecPT = 30.

# RootTupleMaker
process.treeCreator = cms.EDAnalyzer('RootTupleMakerPAT')
process.treeCreator.maxgenparticles = cms.untracked.int32(25)
process.treeCreator.maxgenjets      = cms.untracked.int32(10)
process.treeCreator.maxelectrons    = cms.untracked.int32(10)
process.treeCreator.maxsuperclusters = cms.untracked.int32(10)
process.treeCreator.maxcalojets     = cms.untracked.int32(10)
process.treeCreator.maxmuons        = cms.untracked.int32(10)
process.treeCreator.debug           = cms.untracked.bool(False)
# Overall luminosity normalization  (in pb-1)   
process.treeCreator.luminosity      = cms.untracked.double(100)
process.treeCreator.numEvents       = cms.untracked.int32(10)
process.treeCreator.saveTrigger     = cms.untracked.bool(True)
process.treeCreator.usePDFweight    = cms.untracked.bool(False)
process.treeCreator.PDFSet          = cms.untracked.string("/cteq61.LHgrid")
process.treeCreator.doBeamSpotCorr  = cms.untracked.bool(False)
process.treeCreator.muonLabel       = cms.untracked.InputTag("cleanLayer1Muons");
process.treeCreator.electronLabel   = cms.untracked.InputTag("cleanLayer1Electrons");
process.treeCreator.superClusterLabel = cms.untracked.InputTag("hybridSuperClusters");
process.treeCreator.caloJetLabel    = cms.untracked.InputTag("cleanLayer1Jets");
process.treeCreator.genJetLabel     = cms.untracked.InputTag("iterativeCone5GenJets");
process.treeCreator.electronPt      = cms.untracked.double(30.);
process.treeCreator.electronIso     = cms.untracked.double(0.1);
process.treeCreator.muonPt          = cms.untracked.double(20.);
process.treeCreator.muonIso         = cms.untracked.double(0.05);

# PAT sequence modification
process.patDefaultSequence.remove( process.countLayer1Objects )

#process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path( process.LJFilter*process.patDefaultSequence*process.treeCreator )

# Output module configuration (to enable the PATtuple output, uncomment the lines below)
#process.out = cms.OutputModule("PoolOutputModule",
    #fileName = cms.untracked.string('PATtuple.root'),
    #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    #dropMetaDataForDroppedData = cms.untracked.bool(True),
    #outputCommands = cms.untracked.vstring('drop *')
#)
#process.outpath = cms.EndPath(process.out)
## save PAT Layer 1 output
#from PhysicsTools.PatAlgos.patEventContent_cff import *
#process.out.outputCommands += patEventContent
#process.out.outputCommands += [
    ## GEN
    #'keep recoGenParticles_genParticles_*_*',
    #'keep GenEventInfoProduct_generator_*_*',
    #'keep *_genMetTrue_*_*',
    #'keep *_iterativeCone5GenJets_*_*',
    ## TRIGGER
    #'keep edmTriggerResults_TriggerResults_*_HLT',
    #'keep *_hltTriggerSummaryAOD_*_*',
    ## PAT
    #'keep *_layer1METs*_*_*',
    ## PAT (dropped)
    #'drop *_cleanLayer1Photons_*_*',
    #'drop *_cleanLayer1Taus_*_*', 
    #'drop *_cleanLayer1Hemispheres_*_*' ]
