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

class displacedMuonsPATAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit displacedMuonsPATAnalyzer(const edm::ParameterSet&);
      ~displacedMuonsPATAnalyzer();
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
      edm::EDGetTokenT<edm::View<pat::Muon> > patToken_;
      edm::Handle<edm::View<pat::Muon> > patCollection_;

      //
      // --- Variables used
      //

      // Event info
      Int_t eventId = 0;
      Int_t luminosityBlock = 0;
      Int_t run = 0;

      // Histograms definition

      TH2F* h_PATsegments_DSAsegments;

      // Output definition
      TFile *file_out;

};




////////
//////// -- Constructor
////////
displacedMuonsPATAnalyzer::displacedMuonsPATAnalyzer(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");

   parameters = iConfig;
   
   patToken_  = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));


   h_PATsegments_DSAsegments   = new TH2F("h_PATsegments_DSAsegmentsA", "; muon.numberOfMatches(); nsegments", 15, 0, 15, 15, 0, 15);


}



////////
//////// -- Destructor
////////
displacedMuonsPATAnalyzer::~displacedMuonsPATAnalyzer()
{

}

////////
//////// -- BeginJob
////////
void displacedMuonsPATAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  output_filename = parameters.getParameter<std::string>("nameOfOutput");
  file_out = new TFile(output_filename.c_str(), "RECREATE");

}

////////
//////// -- EndJob
////////
void displacedMuonsPATAnalyzer::endJob()
{
  std::cout << "End Job" << std::endl;
  file_out->cd();
  h_PATsegments_DSAsegments->Write();
  file_out->Close();

}


////////
//////// -- fillDescriptions
////////
void displacedMuonsPATAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


////////
//////// -- Analyze
////////
void displacedMuonsPATAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
     const pat::Muon& muon(patCollection_->at(i));
     std::cout << "Best muon pt: " << muon.pt() << std::endl;
     std::cout << "Best muon inner x: " << muon.bestTrack()->extra()->innerPosition().x() << std::endl;
     if (muon.isGlobalMuon()){
         const reco::Track* dgl = muon.combinedMuon().get();
         //std::cout << "Combined muon pt: " << dgl->pt() << std::endl;
         //std::cout << "Combined muon inner x: " << dgl->innerPosition().x() << std::endl;
     }
     if (muon.isStandAloneMuon()){
       //const reco::Track* dsa = muon.standAloneMuon().get();
       const reco::Track* best = muon.bestTrack();
       //std::cout << "Standalone muon pt: " << dsa->pt() << std::endl;
       unsigned int nbestsegments = 0;
       for (trackingRecHit_iterator hit = best->extra()->recHitsBegin(); hit != best->extra()->recHitsEnd(); ++hit) {
         if (!(*hit)->isValid()) continue;
         DetId id = (*hit)->geographicalId();
         if (id.det() != DetId::Muon) continue;
         if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
           nbestsegments++;
         }
       }

       std::cout << nbestsegments << std::endl;

       h_PATsegments_DSAsegments->Fill(muon.numberOfMatches(), nbestsegments);

     }
   } 
   



}


DEFINE_FWK_MODULE(displacedMuonsPATAnalyzer);
