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

   h_trackerMuons_size               = new TH1F("h_trackerMuons_size", "; Number of Tracker Muons; Events", 10, 0, 10);
   h_trackerMuons_numberOfMatches    = new TH1F("h_trackerMuons_numberOfMatches", "; muon.numberOfMatches(); Events", 14, 0, 14);
   h_trackerMuons_numberOfSegments   = new TH1F("h_trackerMuons_numberOfSegments", "; Number of segments; Events", 30, 0, 30);
   h_trackerMuons_segmentX           = new TH1F("h_trackerMuons_segmentX", "; Segment X position; Events", 50, -200 , 200);
   h_trackerMuons_segmentY           = new TH1F("h_trackerMuons_segmentY", "; Segment Y position; Events", 50, -1 , 1);

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
     }
   }
   h_trackerMuons_size->Fill(nTRK);



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
