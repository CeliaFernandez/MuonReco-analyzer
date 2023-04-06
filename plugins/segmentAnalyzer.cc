#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

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

#include "helper.h"

class segmentAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit segmentAnalyzer(const edm::ParameterSet&);
      ~segmentAnalyzer();
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

      // Muon collection
      edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
      edm::Handle<edm::View<reco::Muon> > muonCollection_;
      edm::EDGetTokenT<edm::View<reco::Track> > trackToken_;
      edm::Handle<edm::View<reco::Track> > trackCollection_;

      // Segment labels
      edm::InputTag theDTRecSegment4DCollectionLabel;
      edm::InputTag theCSCSegmentCollectionLabel;
      edm::InputTag theGEMSegmentCollectionLabel;
      edm::InputTag theRPCHitCollectionLabel;
      edm::InputTag theGEMHitCollectionLabel;
      //edm::InputTag theME0SegmentCollectionLabel;

      // Segment tokens
      edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken;
      edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken;
      edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken;
      edm::EDGetTokenT<RPCRecHitCollection> rpcHitsToken;
      edm::EDGetTokenT<GEMRecHitCollection> gemHitsToken;
      //edm::EDGetTokenT<ME0SegmentCollection> me0SegmentsToken;



      //
      // --- Variables used
      //

      // Event info
      Int_t eventId = 0;
      Int_t luminosityBlock = 0;
      Int_t run = 0;

      // Segments histograms
      // (study the performance of reco segments)
      TH1F* h_total1_number;
      TH1F* h_total2_number;
      TH1F* h_dtSegments_number;
      TH1F* h_cscSegments_number;
      TH1F* h_gemSegments_number;
      TH1F* h_rpcHits_number;
      TH1F* h_gemHits_number;
      TH2F* h_dtSegments_cscSegments_number;

      // Muon histogras histograms
      // (study the performance of tracker muon segments)
      TH1F* h_trackerMuons_size;
      TH1F* h_trackerMuons_numberOfMatches;
      TH1F* h_trackerMuons_numberOfSegments;
      TH1F* h_trackerMuons_segmentX;
      TH1F* h_trackerMuons_segmentY;

      TH1F* h_inner_pt;
      TH1F* h_inner_eta;
      TH1F* h_inner_hits;
      TH1F* h_inner_missingouterhits;
      TH1F* h_inner_missingfraction;
      TH1F* h_inner_chi2;
      TH1F* h_inner_algo;
      TH1F* h_inner_ptError;
      TH1F* h_inner_etaError;
      TH1F* h_inner_phiError;

      TH1F* h_track_pt;
      TH1F* h_track_eta;
      TH1F* h_track_hits;
      TH1F* h_track_missingouterhits;
      TH1F* h_track_missingfraction;
      TH1F* h_track_chi2;
      TH1F* h_track_algo;
      TH1F* h_track_ptError;
      TH1F* h_track_etaError;
      TH1F* h_track_phiError;

      // Output definition
      TFile *file_out;

};




////////
//////// -- Constructor
////////
segmentAnalyzer::segmentAnalyzer(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");

   parameters = iConfig;
   
   muonToken_     = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));
   trackToken_     = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("trackCollection"));

   theDTRecSegment4DCollectionLabel = iConfig.getParameter<edm::InputTag>("DTRecSegment4DCollectionLabel");
   theCSCSegmentCollectionLabel     = iConfig.getParameter<edm::InputTag>("CSCSegmentCollectionLabel");
   theGEMSegmentCollectionLabel     = iConfig.getParameter<edm::InputTag>("GEMSegmentCollectionLabel");
   theRPCHitCollectionLabel         = iConfig.getParameter<edm::InputTag>("RPCHitCollectionLabel");
   theGEMHitCollectionLabel         = iConfig.getParameter<edm::InputTag>("GEMHitCollectionLabel");
   //theME0SegmentCollectionLabel = iConfig.getParameter<edm::InputTag>("ME0SegmentCollectionLabel");

   dtSegmentsToken  = consumes<DTRecSegment4DCollection>(theDTRecSegment4DCollectionLabel);
   cscSegmentsToken = consumes<CSCSegmentCollection>(theCSCSegmentCollectionLabel);
   gemSegmentsToken = consumes<GEMSegmentCollection>(theGEMSegmentCollectionLabel);
   rpcHitsToken     = consumes<RPCRecHitCollection>(theRPCHitCollectionLabel);
   gemHitsToken     = consumes<GEMRecHitCollection>(theGEMHitCollectionLabel);
   //me0SegmentsToken = iC.consumes<ME0SegmentCollection>(theME0SegmentCollectionLabel);

   // Histograms
   h_total1_number        = new TH1F("h_total1_number", "; dtSegments.size() + cscSegments.size() + rpcRecHits.size(); Events", 100, 0, 100);
   h_total2_number        = new TH1F("h_total2_number", "; dtSegments.size() + cscSegments.size() + rpcRecHits.size() + gemSegments.size() + gemRecHits.size(); Events", 100, 0, 100);
   h_dtSegments_number    = new TH1F("h_dtSegments_number", "; dtSegments.size(); Events", 30, 0, 30);
   h_cscSegments_number   = new TH1F("h_cscSegments_number", "; cscSegments.size(); Events", 30, 0, 30);
   h_gemSegments_number   = new TH1F("h_gemSegments_number", "; gemSegments.size(); Events", 30, 0, 30);
   h_rpcHits_number       = new TH1F("h_rpcHits_number", "; rpcRecHits.size(); Events", 30, 0, 30);
   h_gemHits_number       = new TH1F("h_gemHits_number", "; gemRecHits.size(); Events", 30, 0, 30);
   //h_me0Segments_number   = new TH1F("h_me0Segments_number", "; me0Segments.size(); Events", 30, 0, 30);
   h_dtSegments_cscSegments_number   = new TH2F("h_dtSegments_cscSegments_number", "; dtSegments.size(); cscSegments.size()", 30, 0, 30, 30, 0, 30);

   h_trackerMuons_size               = new TH1F("h_trackerMuons_size", "; Number of Tracker Muons; Events", 25, 0, 25);
   h_trackerMuons_numberOfMatches    = new TH1F("h_trackerMuons_numberOfMatches", "; muon.numberOfMatches(); Events", 14, 0, 14);
   h_trackerMuons_numberOfSegments   = new TH1F("h_trackerMuons_numberOfSegments", "; Number of segments; Events", 30, 0, 30);
   h_trackerMuons_segmentX           = new TH1F("h_trackerMuons_segmentX", "; Segment X position; Events", 50, -200 , 200);
   h_trackerMuons_segmentY           = new TH1F("h_trackerMuons_segmentY", "; Segment Y position; Events", 50, -1 , 1);

   // TrackerMuons vs generalTracks
   h_inner_pt  = new TH1F("h_inner_pt", "; Track p_{T}; Events", 80, 0, 80);
   h_inner_eta = new TH1F("h_inner_eta", "; Track #eta; Events", 50, -3, 3);
   h_inner_hits = new TH1F("h_inner_hits", "; Track hits; Events", 50, 0, 50);
   h_inner_missingouterhits = new TH1F("h_inner_missingouterhits", "; Track missing outer hits; Events", 50, 0, 50);
   h_inner_missingfraction = new TH1F("h_inner_missingfraction", "; Track missing hit fraction; Events", 50, 0, 1);
   h_inner_chi2 = new TH1F("h_inner_chi2", "; Track #Chi^2/ndof; Events", 50, 0, 50);
   h_inner_algo = new TH1F("h_inner_algo", "; Track algorithm; Events", 46, 0, 46);
   h_inner_ptError  = new TH1F("h_inner_ptError", "; Track p_{T} error; Events", 50, 0, 50);
   h_inner_etaError  = new TH1F("h_inner_etaError", "; Track #eta error; Events", 50, 0, 50);
   h_inner_phiError  = new TH1F("h_inner_phiError", "; Track #phi error; Events", 50, 0, 50);

   h_track_pt  = new TH1F("h_track_pt", "; Track p_{T}; Events", 80, 0, 80);
   h_track_eta = new TH1F("h_track_eta", "; Track #eta; Events", 50, -3, 3);
   h_track_hits = new TH1F("h_track_hits", "; Track hits; Events", 50, 0, 50);
   h_track_missingouterhits = new TH1F("h_track_missingouterhits", "; Track missing outer hits; Events", 50, 0, 50);
   h_track_missingfraction = new TH1F("h_track_missingfraction", "; Track missing hit fraction; Events", 50, 0, 1);
   h_track_chi2 = new TH1F("h_track_chi2", "; Track #Chi^2/ndof; Events", 50, 0, 50);
   h_track_algo = new TH1F("h_track_algo", "; Track algorithm; Events", 46, 0, 46);
   h_track_ptError  = new TH1F("h_track_ptError", "; Track p_{T} error; Events", 50, 0, 50);
   h_track_etaError  = new TH1F("h_track_etaError", "; Track #eta error; Events", 50, 0, 50);
   h_track_phiError  = new TH1F("h_track_phiError", "; Track #phi error; Events", 50, 0, 50);

}



////////
//////// -- Destructor
////////
segmentAnalyzer::~segmentAnalyzer()
{

}

////////
//////// -- BeginJob
////////
void segmentAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  output_filename = parameters.getParameter<std::string>("nameOfOutput");
  file_out = new TFile(output_filename.c_str(), "RECREATE");

}

////////
//////// -- EndJob
////////
void segmentAnalyzer::endJob()
{
  std::cout << "End Job" << std::endl;
  file_out->cd();

  h_trackerMuons_size->Write();
  h_trackerMuons_numberOfMatches->Write();
  h_trackerMuons_numberOfSegments->Write();
  h_trackerMuons_segmentX->Write();
  h_trackerMuons_segmentY->Write();

  h_total1_number->Write();
  h_total2_number->Write();
  h_dtSegments_number->Write();
  h_cscSegments_number->Write();
  h_gemSegments_number->Write();
  h_rpcHits_number->Write();
  h_gemHits_number->Write();
  h_dtSegments_cscSegments_number->Write();

  h_inner_pt->Write();
  h_inner_eta->Write();
  h_inner_hits->Write();
  h_inner_missingouterhits->Write();
  h_inner_missingfraction->Write();
  h_inner_chi2->Write();
  h_inner_algo->Write();
  h_inner_ptError->Write();
  h_inner_etaError->Write();
  h_inner_phiError->Write();
  h_track_pt->Write();
  h_track_eta->Write();
  h_track_hits->Write();
  h_track_missingfraction->Write();
  h_track_chi2->Write();
  h_track_algo->Write();
  h_track_ptError->Write();
  h_track_etaError->Write();
  h_track_phiError->Write();

  file_out->Close();

}


////////
//////// -- fillDescriptions
////////
void segmentAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


////////
//////// -- Analyze
////////
void segmentAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   //
   // -- Get the collections
   //

   iEvent.getByToken(muonToken_, muonCollection_);
   iEvent.getByToken(trackToken_, trackCollection_);

   edm::Handle<DTRecSegment4DCollection> dtSegments;
   iEvent.getByToken(dtSegmentsToken, dtSegments);

   edm::Handle<CSCSegmentCollection> cscSegments;
   iEvent.getByToken(cscSegmentsToken, cscSegments);

   edm::Handle<GEMSegmentCollection> gemSegments;
   iEvent.getByToken(gemSegmentsToken, gemSegments);

   edm::Handle<RPCRecHitCollection> rpcRecHits;
   iEvent.getByToken(rpcHitsToken, rpcRecHits);

   edm::Handle<GEMRecHitCollection> gemRecHits;
   iEvent.getByToken(gemHitsToken, gemRecHits);


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
   //

   // Fill muon analyiss
   int nTRK;
   int nChambersDT;
   int nChambersCSC;
   int nSegments;
   nTRK = 0;
   for (size_t i = 0; i < muonCollection_->size(); i++) {
     const reco::Muon &muon = (*muonCollection_)[i];
     if (muon.isTrackerMuon()) {

       nTRK++;

       h_trackerMuons_numberOfMatches->Fill(muon.numberOfMatches());
       // Loop over chambers
       nChambersDT = 0;
       nChambersCSC = 0;
       nSegments = 0;
       for (auto& chamber : muon.matches()) {
         if (chamber.detector() == MuonSubdetId::DT)
           nChambersDT++; 
         if (chamber.detector() == MuonSubdetId::CSC)
           nChambersCSC++; 
         for (auto& segment : chamber.segmentMatches) {
           nSegments++;
           h_trackerMuons_segmentX->Fill(segment.x);
           h_trackerMuons_segmentY->Fill(segment.y);
         }
       }
       h_trackerMuons_numberOfSegments->Fill(nSegments);

       if(muon.innerTrack().isNonnull()) {
         h_inner_pt->Fill(muon.innerTrack()->pt());
         h_inner_eta->Fill(muon.innerTrack()->eta());
         h_inner_hits->Fill(muon.innerTrack()->numberOfValidHits());
         h_inner_missingouterhits->Fill(muon.innerTrack()->missingOuterHits());
         h_inner_missingfraction->Fill(muon.innerTrack()->missingOuterHits()/muon.innerTrack()->numberOfValidHits());
         h_inner_chi2->Fill(muon.innerTrack()->normalizedChi2());
         h_inner_algo->Fill(muon.innerTrack()->algo());
         h_inner_ptError->Fill(muon.innerTrack()->ptError());
         h_inner_etaError->Fill(muon.innerTrack()->etaError());
         h_inner_phiError->Fill(muon.innerTrack()->phiError());

       }


     }
   }
   h_trackerMuons_size->Fill(nTRK);

   for (size_t i = 0; i < trackCollection_->size(); i++) {
     const reco::Track &track = trackCollection_->at(i);

     bool isMuon = false; 
     for (size_t j = 0; j < muonCollection_->size(); j++) {
       const reco::Muon &ref = (*muonCollection_)[j];
       if (ref.innerTrack().get() == &track) {
         isMuon = true;
         break;
       }
     }
     if (isMuon)
       continue;
     
     h_track_pt->Fill(track.pt());
     h_track_eta->Fill(track.eta());
     h_track_hits->Fill(track.numberOfValidHits());
     h_track_missingouterhits->Fill(track.missingOuterHits());
     h_track_missingfraction->Fill(track.missingOuterHits()/track.numberOfValidHits());
     h_track_chi2->Fill(track.normalizedChi2());
     h_track_algo->Fill(track.algo());
     h_track_ptError->Fill(track.ptError());
     h_track_etaError->Fill(track.etaError());
     h_track_phiError->Fill(track.phiError());

   }

   // Fill segments histograms
   h_total1_number->Fill(dtSegments->size() + cscSegments->size() + rpcRecHits->size());
   h_total2_number->Fill(dtSegments->size() + cscSegments->size() + rpcRecHits->size() + gemSegments->size() + gemRecHits->size());
   h_dtSegments_number->Fill(dtSegments->size());
   h_cscSegments_number->Fill(cscSegments->size());
   h_gemSegments_number->Fill(gemSegments->size());
   h_rpcHits_number->Fill(rpcRecHits->size());
   h_gemHits_number->Fill(gemRecHits->size());
   h_dtSegments_cscSegments_number->Fill(dtSegments->size(), cscSegments->size());



}


DEFINE_FWK_MODULE(segmentAnalyzer);
