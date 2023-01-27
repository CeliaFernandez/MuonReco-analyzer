#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"

#include "helper.h"

class standardRECOAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit standardRECOAnalyzer(const edm::ParameterSet&);
      ~standardRECOAnalyzer();
      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;
      std::string output_filename;

      //
      // --- Tokens and Handles
      //

      edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
      edm::Handle<edm::View<reco::Muon> > muonCollection_;

      edm::EDGetTokenT<edm::View<reco::GenParticle> > GenToken_;
      edm::Handle<edm::View<reco::GenParticle> > GenCollection_;

      //
      // --- Variables used
      //

      // Event info
      Int_t eventId = 0;
      Int_t luminosityBlock = 0;
      Int_t run = 0;

      // Histograms definition
      TH1F *h_gen_pt;
      TH1F *h_gen_eta;
      TH1F *h_gen_phi;
      TH1F *h_mu_pt;
      TH1F *h_mu_eta;
      TH1F *h_mu_phi;
      TH1F *h_isoR03_sumPt;
      TH1F *h_isoR03_emEt;
      TH1F *h_isoR03_hadEt;
      TH1F *h_isoR03_hoEt;
      TH1F *h_isoR03_nTracks;
      TH1F *h_isoR03_nJets;
      TH1F *h_isoR03_trackerVetoPt;
      TH1F *h_isoR03_emVetoEt;
      TH1F *h_isoR03_hadVetoEt;
      TH1F *h_isoR03_hoVetoEt;
      TH1F *h_pfIsoR03_sumChargedHadronPt;
      TH1F *h_pfIsoR03_sumNeutralHadronEt;
      TH1F *h_pfIsoR03_sumPhotonEt;
      TH1F *h_pfIsoR03_sumNeutralHadronEtHighThreshold;
      TH1F *h_pfIsoR03_sumPhotonEtHighThreshold;
      TH1F *h_pfIsoR03_sumPUPt;
      TH1F *h_glb_normalizedChi2;
      TH1F *h_glb_numberOfValidMuonHits;
      TH1F *h_glb_numberOfMatchedStations;
      TH1F *h_glb_dxy;
      TH1F *h_glb_dz;
      TH1F *h_glb_trackerLayers;
      TH1F *h_glb_pixelLayers;
      TH1F *h_trk_trackerLayers;
      TH1F *h_trk_pixelLayers;
      TH1F *h_trk_dxy;
      TH1F *h_trk_dz;

      // Efficiencies definition
      // reco (wrt gen):
      TEfficiency *e_mu_pt;
      TEfficiency *e_mu_eta;
      TEfficiency *e_mu_phi;
      TEfficiency *e_glb_pt_PFIsoLoose;
      TEfficiency *e_glb_pt_PFIsoMedium;
      TEfficiency *e_glb_pt_PFIsoTight;
      TEfficiency *e_glb_eta_PFIsoLoose;
      TEfficiency *e_glb_eta_PFIsoMedium;
      TEfficiency *e_glb_eta_PFIsoTight;
      TEfficiency *e_glb_pt_TightID;
      TEfficiency *e_glb_eta_TightID;
      TEfficiency *e_promptglb_pt_TightID;
      TEfficiency *e_promptglb_eta_TightID;
      TEfficiency *e_trk_pt_SoftID;
      TEfficiency *e_trk_eta_SoftID;


      // Output definition
      TFile *file_out;

};




////////
//////// -- Constructor
////////
standardRECOAnalyzer::standardRECOAnalyzer(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");

   parameters = iConfig;
   
   muonToken_  = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));
   GenToken_  = consumes<edm::View<reco::GenParticle> >( parameters.getParameter<edm::InputTag>("GenCollection") );

   h_gen_pt = new TH1F("h_gen_pt", ";Gen muon p_{T} (GeV);", 60, 0, 1000);
   h_gen_eta = new TH1F("h_gen_eta", ";Gen muon #eta;", 60, -3, 3);
   h_gen_phi = new TH1F("h_gen_phi", ";Gen muon #phi;", 60, -3.14, 3.14);


   h_mu_pt = new TH1F("h_mu_pt", ";Reco muon p_{T} (GeV);", 60, 0, 1000);
   h_mu_eta = new TH1F("h_mu_eta", ";Reco muon #eta;", 60, -3, 3);
   h_mu_phi = new TH1F("h_mu_phi", ";Reco muon #phi;", 60, -3.14, 3.14);
   h_isoR03_sumPt = new TH1F("h_isoR03_sumPt", ";muon->isolationR03().sumPt;", 20, 0, 20);
   h_isoR03_emEt = new TH1F("h_isoR03_emEt", ";muon->isolationR03().emEt;", 20, 0, 20);
   h_isoR03_hadEt = new TH1F("h_isoR03_hadEt", ";muon->isolationR03().hadEt;", 20, 0, 20);
   h_isoR03_hoEt = new TH1F("h_isoR03_hoEt", ";muon->isolationR03().hoEt;", 20, 0, 20);
   h_isoR03_nTracks = new TH1F("h_isoR03_nTracks", ";muon->isolationR03().nTracks;", 20, 0, 20);
   h_isoR03_nJets = new TH1F("h_isoR03_nJets", ";muon->isolationR03().nJets;", 20, 0, 20);
   h_isoR03_trackerVetoPt = new TH1F("h_isoR03_trackerVetoPt", ";muon->isolationR03().trackerVetoPt;", 40, 0, 100);
   h_isoR03_emVetoEt = new TH1F("h_isoR03_emVetoEt", ";muon->isolationR03().emVetoEt;", 20, 0, 20);
   h_isoR03_hadVetoEt = new TH1F("h_isoR03_hadVetoEt", ";muon->isolationR03().hadVetoEt;", 20, 0, 20);
   h_isoR03_hoVetoEt = new TH1F("h_isoR03_hoVetoEt", ";muon->isolationR03().hoVetoEt;", 20, 0, 20);
   h_pfIsoR03_sumChargedHadronPt = new TH1F("h_pfIsoR03_sumChargedHadronPt", ";muon->pfIsolationR03().sumChargedHadronPt;", 20, 0, 20);
   h_pfIsoR03_sumNeutralHadronEt = new TH1F("h_pfIsoR03_sumNeutralHadronEt", ";muon->pfIsolationR03().sumNeutralHadronEt;", 20, 0, 20);
   h_pfIsoR03_sumPhotonEt = new TH1F("h_pfIsoR03_sumPhotonEt", ";muon->pfIsolationR03().sumPhotonEt;", 20, 0, 20);
   h_pfIsoR03_sumNeutralHadronEtHighThreshold = new TH1F("h_pfIsoR03_sumNeutralHadronEtHighThreshold", ";muon->pfIsolationR03().sumNeutralHadronEtHighThreshold;", 20, 0, 20);
   h_pfIsoR03_sumPhotonEtHighThreshold = new TH1F("h_pfIsoR03_sumPhotonEtHighThreshold", ";muon->pfIsolationR03().sumPhotonEtHighThreshold;", 20, 0, 20);
   h_pfIsoR03_sumPUPt = new TH1F("h_pfIsoR03_sumPUPt", ";muon->pfIsolationR03().sumPUPt;", 20, 0, 20);

   h_glb_normalizedChi2 = new TH1F("h_glb_normalizedChi2", ";GLB muon #chi^{2}/ndf;", 100, 0, 100);
   h_glb_numberOfValidMuonHits = new TH1F("h_glb_numberOfValidMuonHits", ";GLB muon valid hits;", 50, 0, 55);
   h_glb_numberOfMatchedStations = new TH1F("h_glb_numberOfMatchedStations", ";GLB muon matched stations;", 4, 0, 4);
   h_glb_dxy = new TH1F("h_glb_dxy", ";GLB muon |d_{xy}| (cm);", 50, 0, 2);
   h_glb_dz = new TH1F("h_glb_dz", ";GLB muon |d_{z}| (cm);", 50, 0, 5);
   h_glb_trackerLayers = new TH1F("h_glb_trackerLayers", ";GLB muon number of tracker layers;", 25, 0, 25);
   h_glb_pixelLayers = new TH1F("h_glb_pixelLayers", ";GLB muon number of pixel layers;", 10, 0, 10);

   h_trk_trackerLayers = new TH1F("h_trk_trackerLayers", ";TRK muon number of tracker layers;", 25, 0, 25);
   h_trk_pixelLayers = new TH1F("h_trk_pixelLayers", ";TRK muon number of pixel layers;", 10, 0, 10);
   h_trk_dxy = new TH1F("h_trk_dxy", ";TRK muon |d_{xy}| (cm);", 50, 0, 2);
   h_trk_dz = new TH1F("h_trk_dz", ";TRK muon |d_{z}| (cm);", 50, 0, 5);

   e_mu_pt = new TEfficiency("e_mu_pt", ";Generated p_{T} (GeV); Reconstruction efficiency", 60, 0, 1000);
   e_mu_eta = new TEfficiency("e_mu_eta", ";Generated #eta; Reconstruction efficiency", 60, -2.5, 2.5);
   e_mu_phi = new TEfficiency("e_mu_phi", ";Generated #phi; Reconstruction efficiency", 60, -3.14, 3.14);
   e_glb_pt_PFIsoLoose = new TEfficiency("e_glb_pt_PFIsoLoose", ";GLB muon p_{T} (GeV); PFIsoLoose efficiency", 60, 0, 1000);
   e_glb_eta_PFIsoLoose = new TEfficiency("e_glb_eta_PFIsoLoose", ";GLB muon #eta; PFIsoLoose efficiency", 60, -2.5, 2.5);
   e_glb_pt_PFIsoMedium = new TEfficiency("e_glb_pt_PFIsoMedium", ";GLB muon p_{T} (GeV); PFIsoMedium efficiency", 60, 0, 1000);
   e_glb_eta_PFIsoMedium = new TEfficiency("e_glb_eta_PFIsoMedium", ";GLB muon #eta; PFIsoMedium efficiency", 60, -2.5, 2.5);
   e_glb_pt_PFIsoTight = new TEfficiency("e_glb_pt_PFIsoTight", ";GLB muon p_{T} (GeV); PFIsoTight efficiency", 60, 0, 1000);
   e_glb_eta_PFIsoTight = new TEfficiency("e_glb_eta_PFIsoTight", ";GLB muon #eta; PFIsoTight efficiency", 60, -2.5, 2.5);
   e_glb_pt_TightID = new TEfficiency("e_glb_pt_TightID", ";GLB muon p_{T} (GeV); TightID efficiency", 60, 0, 1000);
   e_glb_eta_TightID = new TEfficiency("e_glb_eta_TightID", ";GLB muon #eta; TightID efficiency", 60, -2.5, 2.5);
   e_promptglb_pt_TightID = new TEfficiency("e_promptglb_pt_TightID", ";GLB muon p_{T} (GeV); TightID efficiency", 60, 0, 1000);
   e_promptglb_eta_TightID = new TEfficiency("e_promptglb_eta_TightID", ";GLB muon #eta; TightID efficiency", 60, -2.5, 2.5);
   e_trk_pt_SoftID = new TEfficiency("e_trk_pt_SoftID", ";TRK muon p_{T} (GeV); SoftID efficiency", 60, 0, 1000);
   e_trk_eta_SoftID = new TEfficiency("e_trk_eta_SoftID", ";TRK muon #eta; SoftID efficiency", 60, -2.5, 2.5);

}



////////
//////// -- Destructor
////////
standardRECOAnalyzer::~standardRECOAnalyzer()
{

}

////////
//////// -- BeginJob
////////
void standardRECOAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  output_filename = parameters.getParameter<std::string>("nameOfOutput");
  file_out = new TFile(output_filename.c_str(), "RECREATE");

}

////////
//////// -- EndJob
////////
void standardRECOAnalyzer::endJob()
{
  std::cout << "End Job" << std::endl;
  file_out->cd();
  h_gen_pt->Write();
  h_gen_eta->Write();
  h_gen_phi->Write();
  h_mu_pt->Write();
  h_mu_eta->Write();
  h_mu_phi->Write();
  h_isoR03_sumPt->Write();
  h_isoR03_emEt->Write();
  h_isoR03_hadEt->Write();
  h_isoR03_hoEt->Write();
  h_isoR03_nTracks->Write();
  h_isoR03_nJets->Write();
  h_isoR03_trackerVetoPt->Write();
  h_isoR03_emVetoEt->Write();
  h_isoR03_hadVetoEt->Write();
  h_isoR03_hoVetoEt->Write();
  h_pfIsoR03_sumChargedHadronPt->Write();
  h_pfIsoR03_sumNeutralHadronEt->Write();
  h_pfIsoR03_sumPhotonEt->Write();
  h_pfIsoR03_sumNeutralHadronEtHighThreshold->Write();
  h_pfIsoR03_sumPhotonEtHighThreshold->Write();
  h_pfIsoR03_sumPUPt->Write();
  h_glb_normalizedChi2->Write();
  h_glb_numberOfValidMuonHits->Write();
  h_glb_dxy->Write();
  h_glb_dz->Write();
  h_glb_numberOfMatchedStations->Write();
  h_glb_trackerLayers->Write();
  h_glb_pixelLayers->Write();
  h_trk_trackerLayers->Write();
  h_trk_pixelLayers->Write();
  h_trk_dxy->Write();
  h_trk_dz->Write();
  e_mu_pt->Write();
  e_mu_eta->Write();
  e_mu_phi->Write();
  e_glb_pt_PFIsoLoose->Write();
  e_glb_eta_PFIsoLoose->Write();
  e_glb_pt_PFIsoMedium->Write();
  e_glb_eta_PFIsoMedium->Write();
  e_glb_pt_PFIsoTight->Write();
  e_glb_eta_PFIsoTight->Write();
  e_glb_pt_TightID->Write();
  e_glb_eta_TightID->Write();
  e_promptglb_pt_TightID->Write();
  e_promptglb_eta_TightID->Write();
  e_trk_pt_SoftID->Write();
  e_trk_eta_SoftID->Write();
  file_out->Close();

}


////////
//////// -- fillDescriptions
////////
void standardRECOAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


////////
//////// -- Analyze
////////
void standardRECOAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   iEvent.getByToken(muonToken_, muonCollection_);
   iEvent.getByToken(GenToken_, GenCollection_ );


   //
   // -- Init the variables
   //

   // -> Event info
   eventId = 0;
   luminosityBlock = 0;
   run = 0;

   //
   // -- Pre-analysis
   //

   // -> Event info
   eventId = iEvent.id().event();
   luminosityBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();


   //
   // -- Main analysis
   //

   for (unsigned int i = 0; i < muonCollection_->size(); i++) {

     const reco::Muon& mu(muonCollection_->at(i));

     h_mu_pt->Fill(mu.pt());
     h_mu_eta->Fill(mu.eta());
     h_mu_phi->Fill(mu.phi());


   } 
   
   //
   // -- Efficiency calculator
   // 

   for (const auto& gp : *GenCollection_) {

     if (fabs(gp.pdgId()) != 13)
       continue;

     if (gp.status() != 1)
       continue;

     h_gen_pt->Fill(gp.pt());
     h_gen_eta->Fill(gp.eta());
     h_gen_phi->Fill(gp.phi());

     float delR = 999.;
     reco::Muon bestMu;
     for (const auto& mu : *muonCollection_){
       float _delPhi = deltaPhi(gp.phi(), mu.phi());
       float _delR = sqrt(_delPhi*_delPhi + (mu.eta() - gp.eta())*(mu.eta() - gp.eta()));
       if (_delR < delR) {
         bestMu = mu;
         delR = _delR;
       }
     }

     if (delR < 0.2) {

       // Reco efficiencies
       e_mu_pt->Fill(true, gp.pt());
       e_mu_eta->Fill(true, gp.eta());
       e_mu_phi->Fill(true, gp.phi());

       if ( bestMu.isGlobalMuon() ) {
         e_glb_pt_PFIsoLoose->Fill(bestMu.passed(reco::Muon::PFIsoLoose), bestMu.pt());
         e_glb_pt_PFIsoMedium->Fill(bestMu.passed(reco::Muon::PFIsoMedium), bestMu.pt());
         e_glb_pt_PFIsoTight->Fill(bestMu.passed(reco::Muon::PFIsoTight), bestMu.pt());
         e_glb_eta_PFIsoLoose->Fill(bestMu.passed(reco::Muon::PFIsoLoose), bestMu.eta());
         e_glb_eta_PFIsoMedium->Fill(bestMu.passed(reco::Muon::PFIsoMedium), bestMu.eta());
         e_glb_eta_PFIsoTight->Fill(bestMu.passed(reco::Muon::PFIsoTight), bestMu.eta());

       }

       if ( bestMu.isIsolationValid() ) {
         h_isoR03_sumPt->Fill(bestMu.isolationR03().sumPt);
         h_isoR03_emEt->Fill(bestMu.isolationR03().emEt);
         h_isoR03_hadEt->Fill(bestMu.isolationR03().hadEt);
         h_isoR03_hoEt->Fill(bestMu.isolationR03().hoEt);
         h_isoR03_nTracks->Fill(bestMu.isolationR03().nTracks);
         h_isoR03_nJets->Fill(bestMu.isolationR03().nJets);
         h_isoR03_trackerVetoPt->Fill(bestMu.isolationR03().trackerVetoPt);
         h_isoR03_emVetoEt->Fill(bestMu.isolationR03().emVetoEt);
         h_isoR03_hadVetoEt->Fill(bestMu.isolationR03().hadVetoEt);
         h_isoR03_hoVetoEt->Fill(bestMu.isolationR03().hoVetoEt);
       }

       if ( bestMu.isPFIsolationValid() ) {

         h_pfIsoR03_sumChargedHadronPt->Fill(bestMu.pfIsolationR03().sumChargedHadronPt);
         h_pfIsoR03_sumNeutralHadronEt->Fill(bestMu.pfIsolationR03().sumNeutralHadronEt);
         h_pfIsoR03_sumPhotonEt->Fill(bestMu.pfIsolationR03().sumPhotonEt);
         h_pfIsoR03_sumNeutralHadronEtHighThreshold->Fill(bestMu.pfIsolationR03().sumNeutralHadronEtHighThreshold);
         h_pfIsoR03_sumPhotonEtHighThreshold->Fill(bestMu.pfIsolationR03().sumPhotonEtHighThreshold);
         h_pfIsoR03_sumPUPt->Fill(bestMu.pfIsolationR03().sumPUPt);

       }

       // Tight ID efficiencies
       if (bestMu.isGlobalMuon() && bestMu.isPFMuon()) {
         h_glb_normalizedChi2->Fill(bestMu.globalTrack()->normalizedChi2());
         h_glb_numberOfValidMuonHits->Fill(bestMu.globalTrack()->hitPattern().numberOfValidMuonHits());
         h_glb_numberOfMatchedStations->Fill(bestMu.numberOfMatchedStations());
         h_glb_dxy->Fill(fabs(bestMu.muonBestTrack()->dxy()));
         h_glb_dz->Fill(fabs(bestMu.muonBestTrack()->dz()));
         h_glb_trackerLayers->Fill(bestMu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
         h_glb_pixelLayers->Fill(bestMu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
         e_glb_pt_TightID->Fill(bestMu.passed(reco::Muon::CutBasedIdTight), bestMu.pt());
         e_glb_eta_TightID->Fill(bestMu.passed(reco::Muon::CutBasedIdTight), bestMu.eta());
         if (fabs(bestMu.muonBestTrack()->dxy()) < 0.2 && fabs(bestMu.muonBestTrack()->dz()) < 0.5 ) {
           e_promptglb_pt_TightID->Fill(bestMu.passed(reco::Muon::CutBasedIdTight), bestMu.pt());
           e_promptglb_eta_TightID->Fill(bestMu.passed(reco::Muon::CutBasedIdTight), bestMu.eta());
         }
       }

       // Tight ID efficiencies
       if (bestMu.isTrackerMuon()) {
         h_trk_trackerLayers->Fill(bestMu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
         h_trk_pixelLayers->Fill(bestMu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
         h_trk_dxy->Fill(fabs(bestMu.innerTrack()->dxy()));
         h_trk_dz->Fill(fabs(bestMu.innerTrack()->dz()));
         e_trk_pt_SoftID->Fill(bestMu.passed(reco::Muon::SoftCutBasedId), bestMu.pt());
         e_trk_eta_SoftID->Fill(bestMu.passed(reco::Muon::SoftCutBasedId), bestMu.eta());
       }

     } else {
       e_mu_pt->Fill(false, gp.pt());
       e_mu_eta->Fill(false, gp.eta());
       e_mu_phi->Fill(false, gp.phi());
     }


   }


}


DEFINE_FWK_MODULE(standardRECOAnalyzer);
