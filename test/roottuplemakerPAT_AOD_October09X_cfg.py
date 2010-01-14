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
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring('file:inputFile.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RootTupleMakerPAT_output.root')
)

# load the standard PAT config
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

## Load additional RECO config
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All'
process.load("Configuration.StandardSequences.MagneticField_cff")

# add tcMET and pfMET
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

# fix InputTags for ECAL IsoDeposits (to work with FastSim samples)
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")
process.gamIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
process.gamIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.coreTools import *
restrictInputToAOD(process, ['All'])

# Skim definition
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
process.LJFilterPAT = process.LJFilter.clone()
# Primary skim 
process.LJFilter.muLabel = 'muons'
process.LJFilter.elecLabel = 'gsfElectrons'
process.LJFilter.jetLabel = 'antikt5CaloJets'
process.LJFilter.counteitherleptontype = False
process.LJFilter.muonsMin = -1
process.LJFilter.electronsMin = -1
process.LJFilter.jetsMin = 1
# Secondary skim
process.LJFilterPAT.counteitherleptontype = False
process.LJFilterPAT.muonsMin = -1
process.LJFilterPAT.useElecID = True
process.LJFilterPAT.elecPT = 30.

# RootTupleMaker
process.treeCreator = cms.EDAnalyzer('RootTupleMakerPAT')
process.treeCreator.maxgenparticles = cms.untracked.int32(25)
process.treeCreator.maxgenjets      = cms.untracked.int32(10)
process.treeCreator.maxelectrons    = cms.untracked.int32(10)
process.treeCreator.maxcalojets     = cms.untracked.int32(10)
process.treeCreator.maxmuons        = cms.untracked.int32(10)
process.treeCreator.maxEBsuperclusters= cms.untracked.int32(20)
process.treeCreator.maxEEsuperclusters= cms.untracked.int32(10)
process.treeCreator.debug           = cms.untracked.bool(False)
# overall luminosity normalization  (in pb-1)   
process.treeCreator.luminosity      = cms.untracked.double(100)
process.treeCreator.numEvents       = cms.untracked.int32(10)
process.treeCreator.saveTrigger     = cms.untracked.bool(True)
process.treeCreator.usePDFweight    = cms.untracked.bool(False)
process.treeCreator.PDFSet          = cms.untracked.string("/cteq61.LHgrid")
process.treeCreator.doBeamSpotCorr  = cms.untracked.bool(False)
process.treeCreator.muonLabel       = cms.untracked.InputTag("cleanLayer1Muons");
process.treeCreator.electronLabel   = cms.untracked.InputTag("cleanLayer1Electrons");
process.treeCreator.caloJetLabel    = cms.untracked.InputTag("cleanLayer1Jets");
process.treeCreator.genJetLabel     = cms.untracked.InputTag("iterativeCone5GenJets");
process.treeCreator.ecalEBLabel     = cms.untracked.InputTag("reducedEcalRecHitsEB");
process.treeCreator.ecalEELabel     = cms.untracked.InputTag("reducedEcalRecHitsEE");
process.treeCreator.superClusterEELabel=cms.untracked.InputTag("multi5x5SuperClustersWithPreshower");
process.treeCreator.superClusterEBLabel=cms.untracked.InputTag("hybridSuperClusters");
process.treeCreator.electronPt      = cms.untracked.double(30.);
process.treeCreator.electronIso     = cms.untracked.double(0.1);
process.treeCreator.muonPt          = cms.untracked.double(20.);
process.treeCreator.muonIso         = cms.untracked.double(0.05);

# PAT sequence modification
process.patDefaultSequence.remove( process.countLayer1Objects )

# Path definition
process.p = cms.Path( process.LJFilter*process.patDefaultSequence*process.treeCreator )
#process.p = cms.Path( process.LJFilter*process.patDefaultSequence*process.LJFilterPAT*process.treeCreator )
