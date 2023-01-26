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

      // Efficiencies definition
      // reco (wrt gen):
      TEfficiency *e_mu_pt;
      TEfficiency *e_mu_eta;
      TEfficiency *e_mu_phi;


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

   e_mu_pt = new TEfficiency("e_mu_pt", ";Generated p_{T} (GeV); Efficiency", 60, 0, 1000);
   e_mu_eta = new TEfficiency("e_mu_eta", ";Generated #eta; Efficiency", 60, -3, 3);
   e_mu_phi = new TEfficiency("e_mu_phi", ";Generated #phi; Efficiency", 60, -3.14, 3.14);

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
  e_mu_pt->Write();
  e_mu_eta->Write();
  e_mu_phi->Write();
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

     if ( mu.isIsolationValid() ) {
       h_isoR03_sumPt->Fill(mu.isolationR03().sumPt);
       h_isoR03_emEt->Fill(mu.isolationR03().emEt);
       h_isoR03_hadEt->Fill(mu.isolationR03().hadEt);
       h_isoR03_hoEt->Fill(mu.isolationR03().hoEt);
       h_isoR03_nTracks->Fill(mu.isolationR03().nTracks);
       h_isoR03_nJets->Fill(mu.isolationR03().nJets);
       h_isoR03_trackerVetoPt->Fill(mu.isolationR03().trackerVetoPt);
       h_isoR03_emVetoEt->Fill(mu.isolationR03().emVetoEt);
       h_isoR03_hadVetoEt->Fill(mu.isolationR03().hadVetoEt);
       h_isoR03_hoVetoEt->Fill(mu.isolationR03().hoVetoEt);
     }

     if ( mu.isPFIsolationValid() ) {

       h_pfIsoR03_sumChargedHadronPt->Fill(mu.pfIsolationR03().sumChargedHadronPt);
       h_pfIsoR03_sumNeutralHadronEt->Fill(mu.pfIsolationR03().sumNeutralHadronEt);
       h_pfIsoR03_sumPhotonEt->Fill(mu.pfIsolationR03().sumPhotonEt);
       h_pfIsoR03_sumNeutralHadronEtHighThreshold->Fill(mu.pfIsolationR03().sumNeutralHadronEtHighThreshold);
       h_pfIsoR03_sumPhotonEtHighThreshold->Fill(mu.pfIsolationR03().sumPhotonEtHighThreshold);
       h_pfIsoR03_sumPUPt->Fill(mu.pfIsolationR03().sumPUPt);

     }

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
       e_mu_pt->Fill(true, gp.pt());
       e_mu_eta->Fill(true, gp.eta());
       e_mu_phi->Fill(true, gp.phi());
     } else {
       e_mu_pt->Fill(false, gp.pt());
       e_mu_eta->Fill(false, gp.eta());
       e_mu_phi->Fill(false, gp.phi());
     }


   }


}


DEFINE_FWK_MODULE(standardRECOAnalyzer);
