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

#include "helper.h"

class displacedMuonsRECOAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit displacedMuonsRECOAnalyzer(const edm::ParameterSet&);
      ~displacedMuonsRECOAnalyzer();
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
      edm::EDGetTokenT<edm::View<reco::Muon> > patToken_;
      edm::Handle<edm::View<reco::Muon> > patCollection_;

      //
      // --- Variables used
      //

      // Event info
      Int_t eventId = 0;
      Int_t luminosityBlock = 0;
      Int_t run = 0;

      // Histograms definition

      TH2F* h_MUONsegments_DSAsegments;
      TH2F* h_MUONbestsegments_DSAsegments;
      TH2F* h_MUONmatches_DSAsegments;
      TH2F* h_MUONmatches_DSAchambers;

      // Output definition
      TFile *file_out;

};




////////
//////// -- Constructor
////////
displacedMuonsRECOAnalyzer::displacedMuonsRECOAnalyzer(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");

   parameters = iConfig;
   
   patToken_  = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));


   h_MUONsegments_DSAsegments   = new TH2F("h_MUONsegments_DSAsegments", "; segments from muon.matches(); nsegments", 15, 0, 15, 15, 0, 15);
   h_MUONbestsegments_DSAsegments   = new TH2F("h_MUONbestsegments_DSAsegments", "; segments from muon.bestTrack(); nsegments", 15, 0, 15, 15, 0, 15);
   h_MUONmatches_DSAsegments   = new TH2F("h_MUONmatches_DSAsegments", "; muon.numberOfMatches(); nsegments", 15, 0, 15, 15, 0, 15);
   h_MUONmatches_DSAchambers   = new TH2F("h_MUONmatches_DSAchambers", "; muon.numberOfMatches(); nchambers", 15, 0, 15, 15, 0, 15);


}



////////
//////// -- Destructor
////////
displacedMuonsRECOAnalyzer::~displacedMuonsRECOAnalyzer()
{

}

////////
//////// -- BeginJob
////////
void displacedMuonsRECOAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  output_filename = parameters.getParameter<std::string>("nameOfOutput");
  file_out = new TFile(output_filename.c_str(), "RECREATE");

}

////////
//////// -- EndJob
////////
void displacedMuonsRECOAnalyzer::endJob()
{
  std::cout << "End Job" << std::endl;
  file_out->cd();
  h_MUONsegments_DSAsegments->Write();
  h_MUONbestsegments_DSAsegments->Write();
  h_MUONmatches_DSAsegments->Write();
  h_MUONmatches_DSAchambers->Write();
  file_out->Close();

}


////////
//////// -- fillDescriptions
////////
void displacedMuonsRECOAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


////////
//////// -- Analyze
////////
void displacedMuonsRECOAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   iEvent.getByToken(patToken_, patCollection_);

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
   for (unsigned int i = 0; i < patCollection_->size(); i++) {
     const reco::Muon& muon(patCollection_->at(i));
     if (muon.isStandAloneMuon()){

       const reco::Track* dsa = muon.standAloneMuon().get();
       //std::cout << "Standalone muon pt: " << dsa->pt() << std::endl;
       unsigned int nsegments = 0;
       unsigned int nchambers = 0;
       std::vector<DetId> matchedchambers;
       for (trackingRecHit_iterator hit = dsa->recHitsBegin(); hit != dsa->recHitsEnd(); ++hit) {
         if (!(*hit)->isValid()) continue;
         DetId id = (*hit)->geographicalId();
         if (id.det() != DetId::Muon) continue;
         if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
           nsegments++;
           if (std::find(matchedchambers.begin(), matchedchambers.end(), id) != matchedchambers.end()){
             continue;
           } else {
             nchambers++;
             matchedchambers.push_back(id);
           }
         }
       }

       unsigned int nmuonchambers = 0;
       unsigned int nmuonsegments = 0;
       for (auto& chamber : muon.matches()) {
         if (!(chamber.detector() == MuonSubdetId::DT || chamber.detector() == MuonSubdetId::CSC))
           continue;
         nmuonchambers++;
         for (auto& segment : chamber.segmentMatches) {
           if (!(segment.isMask(reco::MuonSegmentMatch::BelongsToTrackByCleaning))) continue;
           nmuonsegments++;
         }
       }

       const reco::Track *best = muon.bestTrack();
       unsigned int nbestsegments = 0;
       /*
       for (trackingRecHit_iterator hit = best->extra()->recHitsBegin(); hit != best->extra()->recHitsEnd(); ++hit) {
         if (!(*hit)->isValid()) continue;
         DetId id = (*hit)->geographicalId();
         if (id.det() != DetId::Muon) continue;
         if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
           nbestsegments++;
         }
       }
       */

       h_MUONsegments_DSAsegments->Fill(nmuonsegments, nsegments);
       h_MUONmatches_DSAsegments->Fill(muon.numberOfMatches(), nsegments);
       h_MUONmatches_DSAchambers->Fill(muon.numberOfMatches(), nchambers);
       h_MUONbestsegments_DSAsegments->Fill(nbestsegments, nsegments);

     }
   } 
   



}


DEFINE_FWK_MODULE(displacedMuonsRECOAnalyzer);
