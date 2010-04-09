# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

############## IMPORTANT ########################################
# if you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 100
#################################################################

# Output ROOT tree
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RootTupleMakerPAT_output_DATA.root')
)

# Global tag (make sure it always matches with the global tag used to reconstruct input files)
process.GlobalTag.globaltag = 'GR_R_35X_V6::All'

# Events to process
process.maxEvents.input = 100

# Options and Output Report
process.options.wantSummary = True

# Input files
process.source.fileNames = [
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0139/7C27BE8F-B73E-DF11-9780-002354EF3BE2.root'
]

# Turn off MC matching for the process
removeMCMatching(process, ['All'])

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# Get the 7 TeV GeV jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Summer09_7TeV_ReReco332")

# Switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

# Restrict input to AOD
restrictInputToAOD(process, ['All'])

# HEEPify PAT electrons
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
    eleLabel = cms.InputTag("patElectrons"),
    barrelCuts = cms.PSet(heepBarrelCuts),
    endcapCuts = cms.PSet(heepEndcapCuts)
)

# Add 'heepPatElectrons' in the right place and point 'selectedLayer1Electrons' to them
process.patDefaultSequence.replace( process.patElectrons, process.patElectrons*process.heepPatElectrons )
process.selectedPatElectrons.src = cms.InputTag("heepPatElectrons")

# Electron and jet cleaning deltaR parameters
process.cleanPatElectrons.checkOverlaps.muons.deltaR = 0.3
process.cleanPatJets.checkOverlaps.muons.deltaR = 0.5
process.cleanPatJets.checkOverlaps.electrons.deltaR = 0.5

# Skim definition
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
process.LJFilter.muLabel = 'muons'
process.LJFilter.elecLabel = 'gsfElectrons'
process.LJFilter.jetLabel = 'ak5CaloJets'
process.LJFilter.muPT = 10.
process.LJFilter.elecPT = 30.

# RootTupleMaker
process.treeCreator = cms.EDAnalyzer('RootTupleMakerPAT')
process.treeCreator.maxgenparticles = cms.untracked.int32(25)
process.treeCreator.maxgenjets      = cms.untracked.int32(10)
process.treeCreator.maxelectrons    = cms.untracked.int32(10)
process.treeCreator.maxbarrelsuperclusters = cms.untracked.int32(30)
process.treeCreator.maxendcapsuperclusters = cms.untracked.int32(30)
process.treeCreator.maxcalojets     = cms.untracked.int32(10)
process.treeCreator.maxmuons        = cms.untracked.int32(10)
process.treeCreator.debug           = cms.untracked.bool(False)
# overall luminosity normalization  (in pb-1)   
process.treeCreator.luminosity      = cms.untracked.double(100)
process.treeCreator.numEvents       = cms.untracked.int32(10)
process.treeCreator.saveTrigger     = cms.untracked.bool(True)
process.treeCreator.usePDFweight    = cms.untracked.bool(False)
process.treeCreator.PDFSet          = cms.untracked.string("/cteq61.LHgrid")
process.treeCreator.doBeamSpotCorr  = cms.untracked.bool(False)
process.treeCreator.muonLabel       = cms.untracked.InputTag("cleanPatMuons");
process.treeCreator.electronLabel   = cms.untracked.InputTag("cleanPatElectrons");
process.treeCreator.caloJetLabel    = cms.untracked.InputTag("cleanPatJets");
process.treeCreator.genJetLabel     = cms.untracked.InputTag("ak5GenJets");
process.treeCreator.ecalEBLabel     = cms.untracked.InputTag("reducedEcalRecHitsEB");
process.treeCreator.ecalEELabel     = cms.untracked.InputTag("reducedEcalRecHitsEE");
process.treeCreator.superClusterEELabel=cms.untracked.InputTag("multi5x5SuperClustersWithPreshower");
process.treeCreator.superClusterEBLabel=cms.untracked.InputTag("hybridSuperClusters");
process.treeCreator.electronPt      = cms.untracked.double(30.);
process.treeCreator.electronIso     = cms.untracked.double(0.1);
process.treeCreator.muonPt          = cms.untracked.double(10.);
process.treeCreator.muonIso         = cms.untracked.double(0.05);

# Path definition
process.p = cms.Path(
    process.LJFilter*
    process.patDefaultSequence*
    process.treeCreator
)

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

# Schedule definition
process.schedule = cms.Schedule(process.p)
