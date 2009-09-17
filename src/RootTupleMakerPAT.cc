// -*- C++ -*-
//
// Package:    RootTupleMakerPAT
// Class:      RootTupleMakerPAT
// 
/**\class RootTupleMakerPAT RootTupleMakerPAT.cc RootTuple/RootTupleMakerPAT/src/RootTupleMakerPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ellie Lockner
//  PAT version by: Dinko Ferencek
//         Created:  Tue Oct 21 13:56:04 CEST 2008
// $Id: RootTupleMakerPAT.cc,v 1.1 2009/08/21 19:24:07 ferencek Exp $
//
//

#include "Leptoquarks/RootTupleMakerPAT/interface/RootTupleMakerPAT.h"

//namesapecs
using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class RootTupleMakerPAT : public edm::EDAnalyzer {
   public:
      explicit RootTupleMakerPAT(const edm::ParameterSet&);
      ~RootTupleMakerPAT();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      void CreateParticleTree(const edm::Handle<reco::GenParticleCollection> collection);

      void SetTriggers(const edm::Event& iEvent);

      int singleEleRelHLTCounter;
      int muonHLTCounter;

      // LHAPDF stuff
  
      double xfx(const double &x,const double &Q, int fl){
        double f[13], mx = x, mQ = Q;
        evolvepdf_(mx, mQ, f);
        return f[fl+6];
      };

      // ----------member data ---------------------------
      // read from cfg file
      int                  maxgenparticles_;
      int                  maxgenjets_;
      int                  maxelectrons_;
      int                  maxcalojets_;
      int                  maxmuons_;
      bool                 debug_;
      double               luminosity_;
      int                  numEvents_;              
      bool                 saveTrigger_;
      int                  prescaleSingleEleRel_;
      int                  prescaleMuon_;
      bool                 usePDFweight_;
      bool                 doBeamSpotCorr_;
      std::string          PDFset_;
      edm::InputTag        muonLabel_, electronLabel_, caloJetLabel_, genJetLabel_;
      double               electronPt_, electronIso_, muonPt_, muonIso_;

      //Output RootNtuple
      TTree *              m_tree;

      // Event info
      int                  event;
      int                  runnum;
      Float_t              x1;
      Float_t              x2;
      Float_t              Q;
      Float_t              PDFweight[41];

      // Gen Event Quantities  
      float                m_cross_section;
      float                m_auto_cross_section;
      int                  m_processID;
      float                m_filter_eff;
      float                m_pthat;
      float                m_weight;              

      // GenParticles
      Int_t                m_GenParticleCount;
      Float_t              m_GenParticleVX[MAXGENPARTICLES];              
      Float_t              m_GenParticleVY[MAXGENPARTICLES];              
      Float_t              m_GenParticleVZ[MAXGENPARTICLES];              
      Float_t              m_GenParticleP[MAXGENPARTICLES];              
      Float_t              m_GenParticlePt[MAXGENPARTICLES];              
      Float_t              m_GenParticlePx[MAXGENPARTICLES];              
      Float_t              m_GenParticlePy[MAXGENPARTICLES];              
      Float_t              m_GenParticlePz[MAXGENPARTICLES];              
      Float_t              m_GenParticleE[MAXGENPARTICLES];              
      Float_t              m_GenParticleEta[MAXGENPARTICLES];              
      Float_t              m_GenParticlePhi[MAXGENPARTICLES];              
      Int_t                m_GenParticlePdgId[MAXGENPARTICLES];              
      Int_t                m_GenParticleMotherIndex[MAXGENPARTICLES];
      Int_t                m_GenParticleNumDaught[MAXGENPARTICLES];  
      Int_t                m_GenParticleStatus[MAXGENPARTICLES];

      // Trigger
      TString              aNames[MAXHLTBITS];
      char                 aHLTNames[6000];
      Int_t                hltNamesLen;
      Int_t                hltCount;
      bool                 aHLTResults[MAXHLTBITS];

      // Electrons
      Int_t                eleCount;
      Float_t              eleEta[MAXELECTRONS];
      Float_t              elePhi[MAXELECTRONS];
      Float_t              elePt[MAXELECTRONS];
      Float_t              eleEnergy[MAXELECTRONS];
      Int_t                eleCharge[MAXELECTRONS];
      Float_t              eleCaloEnergy[MAXELECTRONS];

      Float_t              eleHoE[MAXELECTRONS];
      Float_t              eleSigmaEE[MAXELECTRONS];
      Float_t              eleDeltaPhiTrkSC[MAXELECTRONS];
      Float_t              eleDeltaEtaTrkSC[MAXELECTRONS];

      Float_t              eleTrkIso[MAXELECTRONS];
      Float_t              eleEcalIso[MAXELECTRONS];
      Float_t              eleHcalIso[MAXELECTRONS];
      Float_t              eleRelIso[MAXELECTRONS];
      Int_t                elePassIso[MAXELECTRONS];
      Int_t                eleOverlaps[MAXELECTRONS];
      Int_t                eleClassif[MAXELECTRONS];
      Float_t              elePassID[MAXELECTRONS];

      // GenJets
      Int_t                genJetCount;
      Float_t              genJetEta[MAXGENJETS];
      Float_t              genJetPhi[MAXGENJETS];
      Float_t              genJetPt[MAXGENJETS];
      Float_t              genJetEnergy[MAXGENJETS];
      Float_t              genJetEMF[MAXGENJETS];
      Float_t              genJetHADF[MAXGENJETS];

      // CaloJets
      Int_t                caloJetCount;
      Float_t              caloJetEta[MAXCALOJETS];
      Float_t              caloJetPhi[MAXCALOJETS];
      Float_t              caloJetEMF[MAXCALOJETS];
      Float_t              caloJetHADF[MAXCALOJETS];
      Float_t              caloJetPt_raw[MAXCALOJETS];
      Float_t              caloJetEnergy_raw[MAXCALOJETS];
      Float_t              caloJetPt[MAXCALOJETS];
      Float_t              caloJetEnergy[MAXCALOJETS];
      Int_t                caloJetOverlaps[MAXCALOJETS];
      Int_t                caloJetPartonFlavour[MAXCALOJETS];
      Float_t              caloJetTrackCountingHighEffBTag[MAXCALOJETS];
      Float_t              caloJetSimpleSecondaryVertexBTag[MAXCALOJETS];
      Float_t              caloJetSoftMuonByPtBTag[MAXCALOJETS];
      Float_t              caloJetBProbabilityBTag[MAXCALOJETS];

      // Muons
      Int_t                muonCount;
      Float_t              muonEta[MAXMUONS];
      Float_t              muonPhi[MAXMUONS];
      Float_t              muonPt[MAXMUONS];
      Float_t              muonEnergy[MAXMUONS];
      Int_t                muonCharge[MAXMUONS];
      Float_t              muonTrkHits[MAXMUONS];
      Float_t              muonTrkD0[MAXMUONS];
      Float_t              muonTrkD0Error[MAXMUONS];
      Float_t              muonTrkDz[MAXMUONS];
      Float_t              muonTrkDzError[MAXMUONS];
      Float_t              muonTrkIso[MAXMUONS];
      Float_t              muonEcalIso[MAXMUONS];
      Float_t              muonHcalIso[MAXMUONS];
      Float_t              muonHOIso[MAXMUONS];
      Float_t              muonRelIso[MAXELECTRONS];
      Int_t                muonPassIso[MAXMUONS];
      Float_t              muonGlobalChi2[MAXMUONS];
      Int_t                muonPassID[MAXMUONS];

      // MET 
      Float_t              genMET;
      Float_t              genMETPhi;
      Float_t              genSumET;
      Float_t              caloMET;
      Float_t              caloMETPhi;
      Float_t              caloSumET;
      Float_t              tcMET;
      Float_t              tcMETPhi;
      Float_t              tcSumET;

      // TFile service
      edm::Service<TFileService> fs;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RootTupleMakerPAT::RootTupleMakerPAT(const edm::ParameterSet& iConfig)

{
  //get parameters from cfg file
  maxgenparticles_   = iConfig.getUntrackedParameter<int>("maxgenparticles",100); 
  maxgenjets_        = iConfig.getUntrackedParameter<int>("maxgenjets",10); 
  maxelectrons_      = iConfig.getUntrackedParameter<int>("maxelectrons",5); 
  maxcalojets_       = iConfig.getUntrackedParameter<int>("maxcalojets",10); 
  maxmuons_          = iConfig.getUntrackedParameter<int>("maxmuons",5); 

  debug_             = iConfig.getUntrackedParameter<bool>("debug",0);
  luminosity_        = iConfig.getUntrackedParameter<double>("luminosity",100); // pb -1
  numEvents_         = iConfig.getUntrackedParameter<int>("numEvents",100); 

  saveTrigger_           = iConfig.getUntrackedParameter<bool>("saveTrigger",1);
  prescaleSingleEleRel_  = iConfig.getUntrackedParameter<int>("prescaleSingleEleRel",30); 
  prescaleMuon_          = iConfig.getUntrackedParameter<int>("prescaleMuon",30);

  usePDFweight_          = iConfig.getUntrackedParameter<bool>("usePDFweight",1);
  PDFset_                = iConfig.getUntrackedParameter<std::string>("PDFSet");
  
  doBeamSpotCorr_        = iConfig.getUntrackedParameter<bool>("doBeamSpotCorr",1);
  
  muonLabel_     = iConfig.getUntrackedParameter<edm::InputTag>("muonLabel",edm::InputTag("cleanLayer1Muons"));
  electronLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("electronLabel",edm::InputTag("cleanLayer1Electrons"));
  caloJetLabel_  = iConfig.getUntrackedParameter<edm::InputTag>("caloJetLabel",edm::InputTag("cleanLayer1Jets"));
  genJetLabel_   = iConfig.getUntrackedParameter<edm::InputTag>("genJetLabel",edm::InputTag("sisCone5GenJets"));
  
  electronPt_    = iConfig.getUntrackedParameter<double>("electronPt",30.);
  electronIso_   = iConfig.getUntrackedParameter<double>("electronIso",0.1);
  muonPt_        = iConfig.getUntrackedParameter<double>("muonPt",20.);
  muonIso_       = iConfig.getUntrackedParameter<double>("muonIso",0.05);

  //Initialize some variables
  singleEleRelHLTCounter=0;
  muonHLTCounter=0;

  event=-999;
  runnum=-999;
  x1=0;
  x2=0;
  Q=0;
  for (int i=0;i<41;i++) PDFweight[i]=0;

  m_cross_section=-999.;
  m_auto_cross_section=-999.;
  m_processID=-999;
  m_filter_eff=-999.;
  m_pthat=-999.;
  m_weight=-999.;              

  genMET=-999.;
  genMETPhi=-999.;
  caloMET=-999.;
  caloMETPhi=-999.;

}


RootTupleMakerPAT::~RootTupleMakerPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void 
RootTupleMakerPAT::beginJob(const edm::EventSetup&)
{
  m_tree = fs->make<TTree>("RootTupleMakerPAT","RootTupleMakerPAT") ;

  m_tree->Branch("event",&event,"event/I");
  m_tree->Branch("run",&runnum,"runnum/I");
  m_tree->Branch("x1",&x1,"x1/F");
  m_tree->Branch("x2",&x2,"x2/F");
  m_tree->Branch("Q",&Q,"Q/F");
  m_tree->Branch("PDFweight",&PDFweight,"PDFWeight[41]/F");

  m_tree->Branch("processID",&m_processID,"processID/I");
  m_tree->Branch("pthat",&m_pthat,"pthat/F");
//   m_tree->Branch("cross_section",&m_cross_section,"cross_section/F");
//   m_tree->Branch("auto_cross_section",&m_auto_cross_section,"auto_cross_section/F");
//   m_tree->Branch("filter_eff",&m_filter_eff,"filter_eff/F"); 
//   m_tree->Branch("weight",&m_weight,"weight/F");
  
  m_tree->Branch("GenParticleCount",&m_GenParticleCount,"GenParticleCount/I");
  m_tree->Branch("GenParticleE",&m_GenParticleE,"GenParticleE[GenParticleCount]/F");
  m_tree->Branch("GenParticleP",&m_GenParticleP,"GenParticleP[GenParticleCount]/F");
  m_tree->Branch("GenParticlePt",&m_GenParticlePt,"GenParticlePt[GenParticleCount]/F");
  m_tree->Branch("GenParticlePx",&m_GenParticlePx,"GenParticlePz[GenParticleCount]/F");
  m_tree->Branch("GenParticlePy",&m_GenParticlePy,"GenParticlePy[GenParticleCount]/F");
  m_tree->Branch("GenParticlePz",&m_GenParticlePz,"GenParticlePz[GenParticleCount]/F");
  m_tree->Branch("GenParticlePdgId",&m_GenParticlePdgId,"GenParticlePdgId[GenParticleCount]/I");
  m_tree->Branch("GenParticleEta",&m_GenParticleEta,"GenParticleEta[GenParticleCount]/F");
  m_tree->Branch("GenParticlePhi",&m_GenParticlePhi,"GenParticlePhi[GenParticleCount]/F");
  m_tree->Branch("GenParticleVX",&m_GenParticleVX,"GenParticleVX[GenParticleCount]/F");
  m_tree->Branch("GenParticleVY",&m_GenParticleVY,"GenParticleVY[GenParticleCount]/F");
  m_tree->Branch("GenParticleVZ",&m_GenParticleVZ,"GenParticleVZ[GenParticleCount]/F");
  m_tree->Branch("GenParticleMotherIndex",&m_GenParticleMotherIndex,"GenParticleMotherIndex[GenParticleCount]/I");
  m_tree->Branch("GenParticleNumDaught",&m_GenParticleNumDaught,"GenParticleNumDaught[GenParticleCount]/I");
  m_tree->Branch("GenParticleStatus",&m_GenParticleStatus,"GenParticleStatus[GenParticleCount]/I");

  m_tree->Branch("hltCount",&hltCount,"hltCount/I");
  m_tree->Branch("hltNamesLen",&hltNamesLen,"hltNamesLen/I");
  m_tree->Branch("HLTNames",&aHLTNames,"HLTNames[hltNamesLen]/C",6000);
  m_tree->Branch("HLTResults",&aHLTResults,"HLTResults[hltCount]/O");

  m_tree->Branch("eleCount",&eleCount,"eleCount/I");
  m_tree->Branch("eleEta",&eleEta,"eleEta[eleCount]/F");
  m_tree->Branch("elePhi",&elePhi,"elePhi[eleCount]/F");
  m_tree->Branch("elePt",&elePt,"elePt[eleCount]/F");
  m_tree->Branch("eleEnergy",&eleEnergy,"eleEnergy[eleCount]/F");
  m_tree->Branch("eleCharge",&eleCharge,"eleCharge[eleCount]/I");
  m_tree->Branch("eleCaloEnergy",&eleCaloEnergy,"eleCaloEnergy[eleCount]/F");
  m_tree->Branch("eleHoE",&eleHoE,"eleHoE[eleCount]/F");
  m_tree->Branch("eleSigmaEE",&eleSigmaEE,"eleSigmaEE[eleCount]/F");
  m_tree->Branch("eleDeltaPhiTrkSC",&eleDeltaPhiTrkSC,"eleDeltaPhiTrkSC[eleCount]/F");
  m_tree->Branch("eleDeltaEtaTrkSC",&eleDeltaEtaTrkSC,"eleDeltaEtaTrkSC[eleCount]/F");
  m_tree->Branch("eleTrkIso",&eleTrkIso,"eleTrkIso[eleCount]/F");
  m_tree->Branch("eleEcalIso",&eleEcalIso,"eleEcalIso[eleCount]/F");
  m_tree->Branch("eleHcalIso",&eleHcalIso,"eleHcalIso[eleCount]/F");
  m_tree->Branch("eleRelIso",&eleRelIso,"eleRelIso[eleCount]/F");
  m_tree->Branch("elePassIso",&elePassIso,"elePassIso[eleCount]/I");
  m_tree->Branch("eleOverlaps",&eleOverlaps,"eleOverlaps[eleCount]/I");
  m_tree->Branch("eleClassif",&eleClassif,"eleClassif[eleCount]/I");
  m_tree->Branch("elePassID",&elePassID,"elePassID[eleCount]/F");

  m_tree->Branch("genJetCount",&genJetCount,"genJetCount/I");
  m_tree->Branch("genJetEta",&genJetEta,"genJetEta[genJetCount]/F");
  m_tree->Branch("genJetPhi",&genJetPhi,"genJetPhi[genJetCount]/F");
  m_tree->Branch("genJetPt",&genJetPt,"genJetPt[genJetCount]/F");
  m_tree->Branch("genJetEnergy",&genJetEnergy,"genJetEnergy[genJetCount]/F");
  m_tree->Branch("genJetEMF",&genJetEMF,"genJetEMF[genJetCount]/F");
  m_tree->Branch("genJetHADF",&genJetHADF,"genJetHADF[genJetCount]/F");

  m_tree->Branch("caloJetCount",&caloJetCount,"caloJetCount/I");
  m_tree->Branch("caloJetEta",&caloJetEta,"caloJetEta[caloJetCount]/F");
  m_tree->Branch("caloJetPhi",&caloJetPhi,"caloJetPhi[caloJetCount]/F");
  m_tree->Branch("caloJetPt",&caloJetPt,"caloJetPt[caloJetCount]/F");
  m_tree->Branch("caloJetEnergy",&caloJetEnergy,"caloJetEnergy[caloJetCount]/F");
  m_tree->Branch("caloJetPt_raw",&caloJetPt_raw,"caloJetPt_raw[caloJetCount]/F");
  m_tree->Branch("caloJetEnergy_raw",&caloJetEnergy_raw,"caloJedtEnergy_raw[caloJetCount]/F");
  m_tree->Branch("caloJetEMF",&caloJetEMF,"caloJetEMF[caloJetCount]/F");
  m_tree->Branch("caloJetHADF",&caloJetHADF,"caloJetHADF[caloJetCount]/F");
  m_tree->Branch("caloJetOverlaps",&caloJetOverlaps,"caloJetOverlaps[caloJetCount]/I");
  m_tree->Branch("caloJetPartonFlavour",&caloJetPartonFlavour,"caloJetPartonFlavour[caloJetCount]/I");
  m_tree->Branch("caloJetTrackCountingHighEffBTag",&caloJetTrackCountingHighEffBTag,"caloJetTrackCountingHighEffBTag[caloJetCount]/F");
  m_tree->Branch("caloJetSimpleSecondaryVertexBTag",&caloJetSimpleSecondaryVertexBTag,"caloJetSimpleSecondaryVertexBTag[caloJetCount]/F");
  m_tree->Branch("caloJetSoftMuonByPtBTag",&caloJetSoftMuonByPtBTag,"caloJetSoftMuonByPtBTag[caloJetCount]/F");
  m_tree->Branch("caloJetBProbabilityBTag",&caloJetBProbabilityBTag,"caloJetBProbabilityBTag[caloJetCount]/F");

  m_tree->Branch("muonCount",&muonCount,"muonCount/I");
  m_tree->Branch("muonEta",&muonEta,"muonEta[muonCount]/F");
  m_tree->Branch("muonPhi",&muonPhi,"muonPhi[muonCount]/F");
  m_tree->Branch("muonPt",&muonPt,"muonPt[muonCount]/F");
  m_tree->Branch("muonEnergy",&muonEnergy,"muonEnergy[muonCount]/F");
  m_tree->Branch("muonCharge",&muonCharge,"muonCharge[muonCount]/I");
  m_tree->Branch("muonTrkHits",&muonTrkHits,"muonTrkHits[muonCount]/F");
  m_tree->Branch("muonTrkD0",&muonTrkD0,"muonTrkD0[muonCount]/F");
  m_tree->Branch("muonTrkD0Error",&muonTrkD0Error,"muonTrkD0Error[muonCount]/F");
  m_tree->Branch("muonTrkDz",&muonTrkDz,"muonTrkDz[muonCount]/F");
  m_tree->Branch("muonTrkDzError",&muonTrkDzError,"muonTrkDzError[muonCount]/F");
  m_tree->Branch("muonTrkIso",&muonTrkIso,"muonTrkIso[muonCount]/F");
  m_tree->Branch("muonEcalIso",&muonEcalIso,"muonEcalIso[muonCount]/F");
  m_tree->Branch("muonHcalIso",&muonHcalIso,"muonHcalIso[muonCount]/F");
  m_tree->Branch("muonHOIso",&muonHOIso,"muonHOIso[muonCount]/F");
  m_tree->Branch("muonRelIso",&muonRelIso,"muonRelIso[muonCount]/F");
  m_tree->Branch("muonPassIso",&muonPassIso,"muonPassIso[muonCount]/I");
  m_tree->Branch("muonGlobalChi2",&muonGlobalChi2,"muonGlobalChi2[muonCount]/F");
  m_tree->Branch("muonPassID",&muonPassID,"muonPassID[muonCount]/I");

  m_tree->Branch("genMET",&genMET,"genMET/F");
  m_tree->Branch("genMETPhi",&genMETPhi,"genMETPhi/F");
  m_tree->Branch("genSumET",&genSumET,"genSumET/F");
  m_tree->Branch("caloMET",&caloMET,"caloMET/F");
  m_tree->Branch("caloMETPhi",&caloMETPhi,"caloMETPhi/F");
  m_tree->Branch("caloSumET",&caloSumET,"caloSumET/F");
  m_tree->Branch("tcMET",&tcMET,"tcMET/F");
  m_tree->Branch("tcMETPhi",&tcMETPhi,"tcMETPhi/F");
  m_tree->Branch("tcSumET",&tcSumET,"tcSumET/F");
}

// ------------ method called for each event  ------------
void
RootTupleMakerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(debug_==true)
    cout << "Analyze " << endl;

  // Fill Event info
  event = iEvent.id().event();
  runnum = iEvent.id().run();

  // Fill gen info
  //using HepMCProduct

  // work in CSA08 RECO //////////////////////////////////
  Handle<HepMCProduct> mcHandle;
  iEvent.getByLabel("source", mcHandle );
  const HepMCProduct* mcCollection = mcHandle.failedToGet () ? 0 : &*mcHandle;
    
  int processID = -99;
  double pthat = -99;
  if (mcCollection) {
      const HepMC::GenEvent *genEvt = mcCollection->GetEvent();
      processID = genEvt->signal_process_id();
  }
  m_processID = processID;

  PdfInfo info; //struct to hold PDF info
  edm::Handle<reco::PdfInfo> pdfHandle;
  iEvent.getByLabel("genEventPdfInfo",pdfHandle);
  const PdfInfo* pdfCollection = pdfHandle.failedToGet () ? 0 : &*pdfHandle;
  if (pdfCollection){
    pthat = pdfCollection->scalePDF;
    x1 = pdfCollection->x1;
    x2 = pdfCollection->x2;
    Q = pdfCollection->scalePDF;
    info.id1 = pdfCollection->id1;
    info.id2 = pdfCollection->id2;
    info.x1 = pdfCollection->x1;
    info.x2 = pdfCollection->x2;
    info.scalePDF = pdfCollection->scalePDF;
  }
  m_pthat = pthat;

  float best_fit=0;
  if (usePDFweight_){
     const char *lhaPDFPath = getenv("LHAPATH");
     string pdfSet(lhaPDFPath);
     string::size_type loc = pdfSet.find(":",0);
     if (loc != string::npos) pdfSet = pdfSet.substr(0,loc);
     pdfSet.append(PDFset_);
     initpdfset_((char *)pdfSet.data(), pdfSet.size());
     // loop over all (error) pdfs
     for (int subpdf=0; subpdf<41; subpdf++){
       initpdf_(subpdf);
       if (subpdf == 0){
         best_fit = xfx(info.x1, info.scalePDF, info.id1)*xfx(info.x2, info.scalePDF, info.id2);
         PDFweight[subpdf]=best_fit;
       }
       else{
          PDFweight[subpdf] =xfx(info.x1, info.scalePDF, info.id1)*xfx(info.x2, info.scalePDF, info.id2)/best_fit;
       }
     }
  }

  // ## work in CSA08 RECO
  //   Handle<double> genEventScale; 
  //   iEvent.getByLabel( "genEventScale", genEventScale );
  //   double pthat = *genEventScale; 
  //   cout << pthat << endl;
  // ##

  // ## don't work in CSA08 RECO
  // double filter_eff = -99.; 
  //   Handle<double> genFilterEff; 
  //   iEvent.getByLabel( "genEventRunInfo", "FilterEfficiency", genFilterEff); 
  //   filter_eff = *genFilterEff;
  //   cout << filter_eff << endl; 
  // ## 

  // ## don't work in CSA08 RECO
  //   double cross_section = -99.; 
  //   Handle<double> genCrossSect; 
  //   iEvent.getByLabel( "genEventRunInfo", "PreCalculatedCrossSection", genCrossSect); 
  //   cross_section = *genCrossSect;
  //   cout << cross_section << endl;
  // ## 

  if(debug_==true)
    cout << "gen event info filled" << endl;

  
  /////////// Electrons
  ///////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(electronLabel_, electrons); 

  // Loop over electrons
  eleCount = 0;
  for( std::vector<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end();++electron ) 
    {
      //exit from loop when you reach the required number of electrons
      if(eleCount > maxelectrons_)
        break;
      
      int overlaps = 0;
      
      const reco::CandidatePtrVector & muons = electron->overlaps("muons");
      for (size_t i = 0; i < muons.size(); ++i) {
        // try to convert to pat::Muons
        const pat::Muon *muon = dynamic_cast<const pat::Muon *>(&*muons[i]);
        if (muon) {
           if (muon->isGood(reco::Muon::GlobalMuonPromptTight) 
               && ((muon->trackIso()+muon->ecalIso()+muon->hcalIso())/muon->pt())<muonIso_
               && muon->pt()>muonPt_) overlaps = 1;
        }
      }

      float trkIso = electron->trackIso();
      float ecalIso = electron->ecalIso();
      float hcalIso = electron->hcalIso();
      float pt = electron->pt();
      float relIso = (trkIso+ecalIso+hcalIso)/pt;
      int passIso = 0;
      if (relIso<electronIso_) passIso = 1;
      int passID = 0;
      /* passID for different electron IDs is assigned bitwise
         bit 0: eidRobustLoose
         bit 1: eidRobustTight
         bit 2: eidLoose
         bit 3: eidTight
         bit 4: eidRobustHighEnergy
      */
      if (electron->electronID("eidRobustLoose")>0) passID = passID | 1<<0;
      if (electron->electronID("eidRobustTight")>0) passID = passID | 1<<1;
      if (electron->electronID("eidLoose")>0) passID = passID | 1<<2;
      if (electron->electronID("eidTight")>0) passID = passID | 1<<3;
      if (electron->electronID("eidRobustHighEnergy")>0) passID = passID | 1<<4;
      
     // Set variables in RootNtuple
      eleEta[eleCount]=electron->eta();
      elePhi[eleCount]=electron->phi();
      elePt[eleCount]=pt;
      eleEnergy[eleCount]=electron->energy();
      eleCharge[eleCount]=electron->charge();
      eleCaloEnergy[eleCount]=electron->caloEnergy();

      //////////ID variables
      eleHoE[eleCount]=electron->hadronicOverEm();
      eleSigmaEE[eleCount]=electron->scSigmaEtaEta();
      eleDeltaPhiTrkSC[eleCount]=electron->deltaPhiSuperClusterTrackAtVtx();
      eleDeltaEtaTrkSC[eleCount]=electron->deltaEtaSuperClusterTrackAtVtx();
      eleClassif[eleCount]=electron->classification();
      elePassID[eleCount]=passID;

      //////////Iso variables
      eleTrkIso[eleCount]=trkIso;
      eleEcalIso[eleCount]=ecalIso;
      eleHcalIso[eleCount]=hcalIso;
      eleRelIso[eleCount]=relIso;
      elePassIso[eleCount]=passIso;
      eleOverlaps[eleCount]=overlaps;

      //go to next electron
      eleCount++;
 
    }


  //////////// Gen Particles
  //////////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel ("genParticles", genParticles);
  CreateParticleTree( genParticles );     
  if(debug_==true)
    cout << "gen particles filled" << endl;

//    reco::GenParticleCollection::const_iterator genParts_it;
//   for (genParts_it = genParticles->begin();
//      genParts_it != genParticles->end(); ++genParts_it)
//        {
//       if (abs(genParts_it->pdgId()) == 13) cout << "Found muon" << endl;
//        }

  //////////// GenJets
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(genJetLabel_, genJets);
  
  genJetCount=0;
  for( GenJetCollection::const_iterator genjet = genJets->begin(); genjet != genJets->end(); genjet++ ) 
    {
      //exit from loop when you reach the required number of electrons
      if(genJetCount > maxgenjets_)
        break;

      float EMF = genjet->emEnergy()/genjet->energy();
      float HADF = genjet->hadEnergy()/genjet->energy();

      genJetPt[genJetCount]=genjet->pt();
      genJetPhi[genJetCount]=genjet->phi();
      genJetEta[genJetCount]=genjet->eta();
      genJetEnergy[genJetCount]=genjet->energy();
      genJetEMF[genJetCount]=EMF;
      genJetHADF[genJetCount]=HADF;
      
      genJetCount++;
    }

  if(debug_==true)
    cout << "genJets filled" << endl;

  ////////////// CaloJets
  ////////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<std::vector<pat::Jet> > caloJets;
  iEvent.getByLabel(caloJetLabel_, caloJets); 

  caloJetCount=0;
  for( std::vector<pat::Jet>::const_iterator calojet = caloJets->begin(); calojet != caloJets->end();++calojet ) 
    {
      //exit from loop when you reach the required number of electrons
      if(caloJetCount > maxcalojets_)
        break;
      
      int overlaps = 0;
      /* overlaps with good electrons (with different electron IDs) and muons are handled bitwise
         bit 0: eidRobustLoose
         bit 1: eidRobustTight
         bit 2: eidLoose
         bit 3: eidTight
         bit 4: eidRobustHighEnergy
         bit 5: GlobalMuonPromptTight
      */
      const reco::CandidatePtrVector & electrons = calojet->overlaps("electrons");
      for (size_t i = 0; i < electrons.size(); ++i) {
        // try to convert to pat::Electron
        const pat::Electron *electron = dynamic_cast<const pat::Electron *>(&*electrons[i]);
        if (electron) {
           if (electron->electronID("eidRobustLoose")>0. 
               && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso_
               && electron->pt()>electronPt_) overlaps = overlaps | 1<<0;
           if (electron->electronID("eidRobustTight")>0. 
               && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso_
               && electron->pt()>electronPt_) overlaps = overlaps | 1<<1;
           if (electron->electronID("eidLoose")>0. 
               && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso_
               && electron->pt()>electronPt_) overlaps = overlaps | 1<<2;
           if (electron->electronID("eidTight")>0. 
               && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso_
               && electron->pt()>electronPt_) overlaps = overlaps | 1<<3;
           if (electron->electronID("eidRobustHighEnergy")>0. 
               && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso_
               && electron->pt()>electronPt_) overlaps = overlaps | 1<<4;
        }
      }
      const reco::CandidatePtrVector & muons = calojet->overlaps("muons");
      for (size_t i = 0; i < muons.size(); ++i) {
        // try to convert to pat::Muons
        const pat::Muon *muon = dynamic_cast<const pat::Muon *>(&*muons[i]);
        if (muon) {
           if (muon->isGood(reco::Muon::GlobalMuonPromptTight) 
               && ((muon->trackIso()+muon->ecalIso()+muon->hcalIso())/muon->pt())<muonIso_
               && muon->pt()>muonPt_) overlaps = overlaps | 1<<5;
        }
      }

      caloJetPt[caloJetCount]=calojet->pt();
      caloJetEnergy[caloJetCount]=calojet->energy();
      caloJetPt_raw[caloJetCount]=calojet->correctedJet("raw").pt();
      caloJetEnergy_raw[caloJetCount]=calojet->correctedJet("raw").energy();
      caloJetPhi[caloJetCount]=calojet->phi();
      caloJetEta[caloJetCount]=calojet->eta();
      caloJetEMF[caloJetCount]=calojet->emEnergyFraction();
      caloJetHADF[caloJetCount]=calojet->energyFractionHadronic();
      caloJetOverlaps[caloJetCount]=overlaps;
      caloJetPartonFlavour[caloJetCount]=calojet->partonFlavour();
      caloJetTrackCountingHighEffBTag[caloJetCount]=calojet->bDiscriminator("trackCountingHighEffBJetTags");
      caloJetSimpleSecondaryVertexBTag[caloJetCount]=calojet->bDiscriminator("simpleSecondaryVertexBJetTag");
      caloJetSoftMuonByPtBTag[caloJetCount]=calojet->bDiscriminator("softMuonByPtBJetTags");
      caloJetBProbabilityBTag[caloJetCount]=calojet->bDiscriminator("jetBProbabilityBJetTags");
      
      caloJetCount++;
    }

  if(debug_==true)
    cout << "CaloJets filled" << endl;

  ////////// MET and GenMET
  //////////////////////////////////////////////////////////////////////////////////////////
  Handle<std::vector<pat::MET> > recoMET;
  iEvent.getByLabel("layer1METs", recoMET);

  const pat::MET& recomet = (*recoMET)[0];
  caloMET = recomet.et();
  caloMETPhi = recomet.phi();
  caloSumET = recomet.sumEt();

  genMET = -99;
  genMETPhi = -99;
  genSumET = -99.;

  const reco::GenMET* genmet = recomet.genMET();
  if (genmet!=NULL) {
     genMET = genmet->pt();
     genMETPhi = genmet->phi();
     genSumET = genmet->sumEt();
  }

  Handle<std::vector<pat::MET> > recoMETtc;
  iEvent.getByLabel("layer1METsTC", recoMETtc);

  const pat::MET& recomettc = (*recoMETtc)[0];
  tcMET = recomettc.et();
  tcMETPhi = recomettc.phi();
  tcSumET = recomettc.sumEt();

  if(debug_==true)
    cout << "MET filled" << endl;

  /////////////////// HLT info
  ///////////////////////////////////////////////////////////////////////////////////////////////

  if(saveTrigger_==true)  
    {
      SetTriggers(iEvent);
      if(debug_==true)
        cout << "HLT bits filled " << endl;
    }
  else 
    {
      // reset variables in ntuple
      for(unsigned int iHLT=0; iHLT<MAXHLTBITS; ++iHLT) {
        aHLTResults[iHLT] = false;
      }

      hltCount=-999;
      hltNamesLen=-999;
      strcpy(aHLTNames,"");

      if(debug_==true)
        cout << "HLT bits not filled" << endl;
    }

  ///////////////// Muons + BeamSpot
  /////////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(muonLabel_, muons);
  
  Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot", beamSpot );  

  muonCount = 0;  
  for( std::vector<pat::Muon>::const_iterator muon = muons->begin(); muon != muons->end();++muon )
    {
      //exit from loop when you reach the required number of muons
      if(muonCount > maxmuons_)
        break;

      if(!muon->isGlobalMuon())
        continue;
        
      float trkd0 = muon->track()->d0();
      if (doBeamSpotCorr_ && beamSpot.isValid()) trkd0 = -(muon->track()->dxy( beamSpot->position()));
      else if (doBeamSpotCorr_ && !beamSpot.isValid()) cout << "No beam spot available!" << endl;
      float trkIso = muon->trackIso();
      float ecalIso = muon->ecalIso();
      float hcalIso = muon->hcalIso();
      float pt = muon->pt();
      float relIso = (trkIso+ecalIso+hcalIso)/pt;
      int passIso = 0;
      if (relIso<muonIso_) passIso = 1;
      int passID = 0;
      if (muon->isGood(reco::Muon::GlobalMuonPromptTight)) passID = 1;
      
      muonEta[muonCount]        = muon->eta();
      muonPhi[muonCount]        = muon->phi();
      muonPt[muonCount]         = pt;
      muonEnergy[muonCount]     = muon->energy();
      muonCharge[muonCount]     = muon->charge();
      muonTrkHits[muonCount]    = muon->track()->numberOfValidHits();
      muonTrkD0[muonCount]      = trkd0;
      muonTrkD0Error[muonCount] = muon->track()->d0Error();
      muonTrkDz[muonCount]      = muon->track()->dz();
      muonTrkDzError[muonCount] = muon->track()->dzError();
      muonGlobalChi2[muonCount] = muon->track()->normalizedChi2();
      muonTrkIso[muonCount]     = trkIso;
      muonEcalIso[muonCount]    = ecalIso;
      muonHcalIso[muonCount]    = hcalIso;
      muonHOIso[muonCount]      = muon->isolationR03().hoEt;
      muonRelIso[muonCount]     = relIso;
      muonPassIso[muonCount]    = passIso;
      muonPassID[muonCount]     = passID;

      muonCount++;
    } 

  //////////////////////////////////////////////////////////////////////////////////////////////

  //  Fill tree for each event
  //  ********************************************************
  //  ********************************************************
  if(debug_==true)
    cout << "About to fill tree" << endl;
  
  m_tree->Fill();  
}

////====================================================================
void RootTupleMakerPAT::CreateParticleTree(edm::Handle<reco::GenParticleCollection> collection)
{
////====================================================================
  m_GenParticleCount = 0;
  reco::GenParticleCollection::const_iterator cand;
  for(cand=collection->begin(); cand!=collection->end(); cand++) 
  {

    if(m_GenParticleCount  > maxgenparticles_)
        break;

    m_GenParticleP[m_GenParticleCount] = cand->p();
    m_GenParticlePx[m_GenParticleCount] = cand->px();
    m_GenParticlePy[m_GenParticleCount] = cand->py();
    m_GenParticlePz[m_GenParticleCount] = cand->pz();
    m_GenParticlePt[m_GenParticleCount] = cand->pt();
    m_GenParticleEta[m_GenParticleCount] = cand->eta();
    m_GenParticlePhi[m_GenParticleCount] = cand->phi();
    m_GenParticleE[m_GenParticleCount] = cand->energy();
    m_GenParticlePdgId[m_GenParticleCount] = cand->pdgId();
    m_GenParticleNumDaught[m_GenParticleCount] = cand->numberOfDaughters();
    m_GenParticleStatus[m_GenParticleCount] = cand->status();
   
    m_GenParticleVX[m_GenParticleCount] = cand->vx();
    m_GenParticleVY[m_GenParticleCount] = cand->vy();
    m_GenParticleVZ[m_GenParticleCount] = cand->vz();
 
    m_GenParticleMotherIndex[m_GenParticleCount] = -1;
    int idx=0;
    reco::GenParticleCollection::const_iterator candIter;
    for(candIter = collection->begin(); candIter!=collection->end(); candIter++) {
      if(&(*candIter)==cand->mother()) {
       m_GenParticleMotherIndex[m_GenParticleCount] = idx;
       break;
      }
      idx++;
    }
  m_GenParticleCount++;
  }    
}


////====================================================================
void RootTupleMakerPAT::SetTriggers(const edm::Event& iEvent)
{
  
////====================================================================
  // reset variables in ntuple
  for(unsigned int iHLT=0; iHLT<MAXHLTBITS; ++iHLT) {
    aHLTResults[iHLT] = false;
  }

  strcpy(aHLTNames,"");
  hltNamesLen = 0;

  edm::InputTag tag("TriggerResults::HLT");
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  iEvent.getByLabel(tag,hltTriggerResults);


  if(!hltTriggerResults.isValid()) {
    std::cout << "invalid handle for HLT TriggerResults" << std::endl;
  } else {
    
    edm::TriggerNames HLTNames;
    HLTNames.init(*hltTriggerResults);
    
    std::string tempnames;
    
    hltCount = hltTriggerResults->size();

    for(int i = 0 ; i < hltCount ; i++) {
      //cout << "HLTTrigger: " << HLTNames.triggerName(i).c_str() << "  " << hltTriggerResults->accept(i) << endl;
      tempnames += HLTNames.triggerName(i) + ":";
      aHLTResults[i] = hltTriggerResults->accept(i);
    }
    
    hltNamesLen = tempnames.length();
    strcpy(aHLTNames,tempnames.c_str());
  
  }


}


// ------------ method called once each job just after ending the event loop  ------------
void 
RootTupleMakerPAT::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(RootTupleMakerPAT);
