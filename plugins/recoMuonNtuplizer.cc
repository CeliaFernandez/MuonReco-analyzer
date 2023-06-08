#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"



class recoMuonNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit recoMuonNtuplizer(const edm::ParameterSet&);
      ~recoMuonNtuplizer();

      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;

      //
      // --- Tokens and Handles
      //

      // muons (reco::Muon // pat::Muon)
      edm::EDGetTokenT<edm::View<reco::Muon> > muToken;
      edm::Handle<edm::View<reco::Muon> > muons;

      edm::EDGetTokenT<edm::View<reco::Track> > trackToken;
      edm::Handle<edm::View<reco::Track> > tracks;

      //
      // --- Variables
      //

      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      Int_t ntk = 0;
      Float_t tk_pt[500] = {0.};

      // Muons
      Int_t nmu = 0;
      Float_t mu_pt[200] = {0.};
      Float_t mu_eta[200] = {0.};
      Float_t mu_phi[200] = {0.};
      Int_t mu_isStandAloneMuon[200] = {0};
      Int_t mu_isGlobalMuon[200] = {0};
      Int_t mu_isTrackerMuon[200] = {0};
      Int_t mu_isCaloMuon[200] = {0};
      Int_t mu_isPFMuon[200] = {0};
      Int_t mu_isRPCMuon[200] = {0};
      Int_t mu_isGEMMuon[200] = {0};
      Int_t mu_isME0Muon[200] = {0};
      Int_t mu_isEnergyValid[200] = {0};
      Int_t mu_isQualityValid[200] = {0};

      // Segment info
      Int_t mu_isMatchesValid[200] = {0};
      Int_t mu_numberOfChambers[200] = {0};
      Int_t mu_numberOfChambersCSCorDT[200] = {0};
      Int_t mu_numberOfMatches[200] = {0};
      Int_t mu_numberOfMatchedStations[200] = {0};
      Int_t mu_expectedNnumberOfMatchedStations[200] = {0};
      Int_t mu_numberOfMatchedRPCLayers[200] = {0};
      Int_t mu_nDigisInStation[200] = {0};
      Int_t mu_numberOfShowers[200] = {0};
      Int_t mu_numberOfSegments[200] = {0};


      // Time
      Int_t mu_isTimeValid[200] = {0};
      Float_t mu_time_nDof[200] = {0.};
      Float_t mu_time_timeAtIpInOut[200] = {0.};
      Float_t mu_time_timeAtIpInOutErr[200] = {0.};
      Float_t mu_time_timeAtIpOutIn[200] = {0.};
      Float_t mu_time_timeAtIpOutInErr[200] = {0.};
      Float_t mu_rpctime_nDof[200] = {0.};
      Float_t mu_rpctime_timeAtIpInOut[200] = {0.};
      Float_t mu_rpctime_timeAtIpInOutErr[200] = {0.};
      Float_t mu_rpctime_timeAtIpOutIn[200] = {0.};
      Float_t mu_rpctime_timeAtIpOutInErr[200] = {0.};
   
      // ID
      Int_t mu_selector[36][200];
      std::vector<std::string> Selectors_;
      bool passSelector[36] = {false};
 

      //
      // --- Output
      //
      std::string output_filename;
      TH1F *counts;
      TFile *file_out;
      TTree *tree_out;

};

// Constructor
recoMuonNtuplizer::recoMuonNtuplizer(const edm::ParameterSet& iConfig) {

   usesResource("TFileService");

   parameters = iConfig;

   counts = new TH1F("counts", "", 1, 0, 1);

   muToken = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));
   trackToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("trackCollection"));

}


// Destructor
recoMuonNtuplizer::~recoMuonNtuplizer() {
}


// beginJob (Before first event)
void recoMuonNtuplizer::beginJob() {

   std::cout << "Begin Job" << std::endl;

   // IDs
   Selectors_.push_back("CutBasedIdLoose");
   Selectors_.push_back("CutBasedIdMedium");
   Selectors_.push_back("CutBasedIdMediumPrompt");
   Selectors_.push_back("CutBasedIdTight");
   Selectors_.push_back("CutBasedIdGlobalHighPt");
   Selectors_.push_back("CutBasedIdTrkHighPt");
   Selectors_.push_back("PFIsoVeryLoose");
   Selectors_.push_back("PFIsoLoose");
   Selectors_.push_back("PFIsoMedium");
   Selectors_.push_back("PFIsoTight");
   Selectors_.push_back("PFIsoVeryTight");
   Selectors_.push_back("TkIsoLoose");
   Selectors_.push_back("TkIsoTight");
   Selectors_.push_back("SoftCutBasedId");
   Selectors_.push_back("SoftMvaId");
   Selectors_.push_back("MvaLoose");
   Selectors_.push_back("MvaMedium");
   Selectors_.push_back("MvaTight");
   Selectors_.push_back("MiniIsoLoose");
   Selectors_.push_back("MiniIsoMedium");
   Selectors_.push_back("MiniIsoTight");
   Selectors_.push_back("MiniIsoVeryTight");
   Selectors_.push_back("TriggerIdLoose");
   Selectors_.push_back("InTimeMuon");
   Selectors_.push_back("PFIsoVeryVeryTight");
   Selectors_.push_back("MultiIsoLoose");
   Selectors_.push_back("MultiIsoMedium");
   Selectors_.push_back("PuppiIsoLoose");
   Selectors_.push_back("PuppiIsoMedium");
   Selectors_.push_back("PuppiIsoTight");
   Selectors_.push_back("MvaVTight");
   Selectors_.push_back("MvaVVTight");
   Selectors_.push_back("LowPtMvaLoose");
   Selectors_.push_back("LowPtMvaMedium");
   Selectors_.push_back("MvaIDwpMedium");
   Selectors_.push_back("MvaIDwpTight");

   // Init the file and the TTree
   output_filename = parameters.getParameter<std::string>("nameOfOutput");
   file_out = new TFile(output_filename.c_str(), "RECREATE");
   tree_out = new TTree("Muons", "Muons");

   // TTree branches
   tree_out->Branch("RecoTrack_number", &ntk, "ntk/I");
   tree_out->Branch("RecoTrack_pt", tk_pt, "tk_pt[ntk]/F");

   tree_out->Branch("RecoMuon_number", &nmu, "nmu/I");
   tree_out->Branch("RecoMuon_pt", mu_pt, "mu_pt[nmu]/F");
   tree_out->Branch("RecoMuon_eta", mu_eta, "mu_eta[nmu]/F");
   tree_out->Branch("RecoMuon_phi", mu_phi, "mu_phi[nmu]/F");
   tree_out->Branch("RecoMuon_isStandAloneMuon", mu_isStandAloneMuon, "mu_isStandAloneMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isGlobalMuon", mu_isGlobalMuon, "mu_isGlobalMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isTrackerMuon", mu_isTrackerMuon, "mu_isTrackerMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isCaloMuon", mu_isCaloMuon, "mu_isCaloMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isPFMuon", mu_isPFMuon, "mu_isPFMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isRPCMuon", mu_isRPCMuon, "mu_isRPCMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isGEMMuon", mu_isGEMMuon, "mu_isGEMMuon[nmu]/I");
   tree_out->Branch("RecoMuon_isME0Muon", mu_isME0Muon, "mu_isME0Muon[nmu]/I");
   tree_out->Branch("RecoMuon_isEnergyValid", mu_isEnergyValid, "mu_isEnergyValid[nmu]/I");
   tree_out->Branch("RecoMuon_isQualityValid", mu_isQualityValid, "mu_isQualityValid[nmu]/I");

   tree_out->Branch("RecoMuon_isMatchesValid", mu_isMatchesValid, "mu_isMatchesValid[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfChambers", mu_numberOfChambers, "mu_numberOfChambers[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfChambersCSCorDT", mu_numberOfChambersCSCorDT, "mu_numberOfChambersCSCorDT[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfMatches", mu_numberOfMatches, "mu_numberOfMatches[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfMatchedStations", mu_numberOfMatchedStations, "mu_numberOfMatchedStations[nmu]/I");
   tree_out->Branch("RecoMuon_expectedNnumberOfMatchedStations", mu_expectedNnumberOfMatchedStations, "mu_expectedNnumberOfMatchedStations[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfMatchedRPCLayers", mu_numberOfMatchedRPCLayers, "mu_numberOfMatchedRPCLayers[nmu]/I");
   tree_out->Branch("RecoMuon_nDigisInStation", mu_nDigisInStation, "mu_nDigisInStation[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfShowers", mu_numberOfShowers, "mu_numberOfShowers[nmu]/I");
   tree_out->Branch("RecoMuon_numberOfSegments", mu_numberOfSegments, "mu_numberOfSegments[nmu]/I");

   tree_out->Branch("RecoMuon_isTimeValid", mu_isTimeValid, "mu_isTimeValid[nmu]/I");
   tree_out->Branch("RecoMuon_time_nDof", mu_time_nDof, "mu_time_nDof[nmu]/F");
   tree_out->Branch("RecoMuon_time_timeAtIpInOut", mu_time_timeAtIpInOut, "mu_time_timeAtIpInOut[nmu]/F");
   tree_out->Branch("RecoMuon_time_timeAtIpInOutErr", mu_time_timeAtIpInOutErr, "mu_time_timeAtIpInOutErr[nmu]/F");
   tree_out->Branch("RecoMuon_time_timeAtIpOutIn", mu_time_timeAtIpOutIn, "mu_time_timeAtIpOutIn[nmu]/F");
   tree_out->Branch("RecoMuon_time_timeAtIpOutInErr", mu_time_timeAtIpOutInErr, "mu_time_timeAtIpOutInErr[nmu]/F");
   tree_out->Branch("RecoMuon_rpctime_nDof", mu_rpctime_nDof, "mu_rpctime_nDof[nmu]/F");
   tree_out->Branch("RecoMuon_rpctime_timeAtIpInOut", mu_rpctime_timeAtIpInOut, "mu_rpctime_timeAtIpInOut[nmu]/F");
   tree_out->Branch("RecoMuon_rpctime_timeAtIpInOutErr", mu_rpctime_timeAtIpInOutErr, "mu_rpctime_timeAtIpInOutErr[nmu]/F");
   tree_out->Branch("RecoMuon_rpctime_timeAtIpOutIn", mu_rpctime_timeAtIpOutIn, "mu_rpctime_timeAtIpOutIn[nmu]/F");
   tree_out->Branch("RecoMuon_rpctime_timeAtIpOutInErr", mu_rpctime_timeAtIpOutInErr, "mu_rpctime_timeAtIpOutInErr[nmu]/F");

   // ID
   std::string label = "RecoMuon_";
   std::string olabel = "[nmu]/I";
   for (unsigned int id = 0; id < Selectors_.size(); id++) {
   //  tree_out->Branch(TString(Selectors_[id]), &passSelector[id]);
     tree_out->Branch(label + TString(Selectors_[id]), &mu_selector[id], label + TString(Selectors_[id])+olabel);
   }

}

// endJob (After event loop has finished)
void recoMuonNtuplizer::endJob()
{

    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();

}


// fillDescriptions
void recoMuonNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// Analyze (per event)
void recoMuonNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   iEvent.getByToken(muToken, muons);
   iEvent.getByToken(trackToken, tracks);

   // Count number of events read
   counts->Fill(0);

   // Tracks
   ntk = 0;
   for (unsigned int i = 0; i < tracks->size(); i++) {
     const reco::Track& track(tracks->at(i));

     tk_pt[ntk] = track.pt();

     ntk++;
   }

   // Muons
   nmu = 0;
   for (unsigned int i = 0; i < muons->size(); i++) {
     const reco::Muon& muon(muons->at(i));

     mu_pt[nmu] = muon.pt();
     mu_eta[nmu] = muon.eta();
     mu_phi[nmu] = muon.phi();
     mu_isStandAloneMuon[nmu] = muon.isStandAloneMuon();
     mu_isGlobalMuon[nmu] = muon.isGlobalMuon();
     mu_isTrackerMuon[nmu] = muon.isTrackerMuon();
     mu_isCaloMuon[nmu] = muon.isCaloMuon();
     mu_isPFMuon[nmu] = muon.isPFMuon();
     mu_isRPCMuon[nmu] = muon.isRPCMuon();
     mu_isGEMMuon[nmu] = muon.isGEMMuon();
     mu_isME0Muon[nmu] = muon.isME0Muon();
     mu_isEnergyValid[nmu] = muon.isEnergyValid();
     mu_isQualityValid[nmu] = muon.isQualityValid();

     mu_isMatchesValid[nmu] = muon.isMatchesValid();
     mu_numberOfChambers[nmu] = muon.numberOfChambers();
     mu_numberOfChambersCSCorDT[nmu] = muon.numberOfChambersCSCorDT();
     mu_numberOfMatches[nmu] = muon.numberOfMatches();
     mu_numberOfMatchedStations[nmu] = muon.numberOfMatchedStations();
     mu_expectedNnumberOfMatchedStations[nmu] = muon.expectedNnumberOfMatchedStations();
     mu_numberOfMatchedRPCLayers[nmu] = muon.numberOfMatchedRPCLayers();
     //mu_nDigisInStation[nmu] = muon.nDigisInStation();
     mu_numberOfShowers[nmu] = muon.numberOfShowers();
     //mu_numberOfSegments[nmu] = muon.numberOfSegments();


     
     mu_isTimeValid[nmu] = muon.isTimeValid();
     mu_time_nDof[nmu] = muon.time().nDof;
     mu_time_timeAtIpInOut[nmu] = muon.time().timeAtIpInOut;
     mu_time_timeAtIpInOutErr[nmu] = muon.time().timeAtIpInOutErr;
     mu_time_timeAtIpOutIn[nmu] = muon.time().timeAtIpOutIn;
     mu_time_timeAtIpOutInErr[nmu] = muon.time().timeAtIpOutInErr;
     mu_rpctime_nDof[nmu] = muon.rpcTime().nDof;
     mu_rpctime_timeAtIpInOut[nmu] = muon.rpcTime().timeAtIpInOut;
     mu_rpctime_timeAtIpInOutErr[nmu] = muon.rpcTime().timeAtIpInOutErr;
     mu_rpctime_timeAtIpOutIn[nmu] = muon.rpcTime().timeAtIpOutIn;
     mu_rpctime_timeAtIpOutInErr[nmu] = muon.rpcTime().timeAtIpOutInErr;

     for (unsigned int id = 0; id < Selectors_.size(); id++) {
        //passSelector[id] = muon.passed(muon::selectorFromString(Selectors_.at(id)));
        mu_selector[id][nmu] = muon.passed(muon::selectorFromString(Selectors_.at(id)));
     }


     nmu++;
   }


   // -> Fill tree
   tree_out->Fill();

}

DEFINE_FWK_MODULE(recoMuonNtuplizer);
