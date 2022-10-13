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
      //edm::EDGetTokenT<edm::View<pat::Muon> > patToken_;
      //edm::Handle<edm::View<pat::Muon> > patCollection_;

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

      // Histograms definition

      TH1F* h_total1_number;
      TH1F* h_total2_number;
      TH1F* h_dtSegments_number;
      TH1F* h_cscSegments_number;
      TH1F* h_gemSegments_number;
      TH1F* h_rpcHits_number;
      TH1F* h_gemHits_number;
      TH2F* h_dtSegments_cscSegments_number;

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

   //gemSegmentsToken = iC.consumes<GEMSegmentCollection>(theGEMSegmentCollectionLabel);
   //me0SegmentsToken = iC.consumes<ME0SegmentCollection>(theME0SegmentCollectionLabel);

   h_total1_number        = new TH1F("h_total1_number", "; dtSegments.size() + cscSegments.size() + rpcRecHits.size(); Events", 100, 0, 100);
   h_total2_number        = new TH1F("h_total2_number", "; dtSegments.size() + cscSegments.size() + rpcRecHits.size() + gemSegments.size() + gemRecHits.size(); Events", 100, 0, 100);
   h_dtSegments_number    = new TH1F("h_dtSegments_number", "; dtSegments.size(); Events", 30, 0, 30);
   h_cscSegments_number   = new TH1F("h_cscSegments_number", "; cscSegments.size(); Events", 30, 0, 30);
   h_gemSegments_number   = new TH1F("h_gemSegments_number", "; gemSegments.size(); Events", 30, 0, 30);
   h_rpcHits_number       = new TH1F("h_rpcHits_number", "; rpcRecHits.size(); Events", 30, 0, 30);
   h_gemHits_number       = new TH1F("h_gemHits_number", "; gemRecHits.size(); Events", 30, 0, 30);
   //h_me0Segments_number   = new TH1F("h_me0Segments_number", "; me0Segments.size(); Events", 30, 0, 30);
   h_dtSegments_cscSegments_number   = new TH2F("h_dtSegments_cscSegments_number", "; dtSegments.size(); cscSegments.size()", 30, 0, 30, 30, 0, 30);


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
