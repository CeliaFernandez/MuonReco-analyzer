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

class standaloneDuplicates : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit standaloneDuplicates(const edm::ParameterSet&);
      ~standaloneDuplicates();
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

      // Gen collection
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
      edm::Handle<edm::View<reco::GenParticle> > genCollection_;

      //
      // --- Variables used
      //

      // Event info
      Int_t eventId = 0;
      Int_t luminosityBlock = 0;
      Int_t run = 0;

      // Histograms definition

      TH1F* h_nSTA;
      TH1F* h_STA_pt;
      TH1F* h_STA_eta;
      TH1F* h_STA_phi;

      TH1F* h_nPairs;
      TH1F* h_Nextra;
      TH1F* h_Nshared;
      TH1F* h_nGLBsInPair;
      TH1F* h_nCommonSegments;
      TH2F* h_nCommonSegments_nCommonHits;
      TH1F* h_nCommonHits;

      TH2F* h_ptGLB_ptSTA;           
      TH2F* h_etaGLB_etaSTA;           
      TH2F* h_phiGLB_phiSTA;           
      TH2F* h_nDTGLB_nDTSTA;           
      TH2F* h_nCSCGLB_nCSCSTA;           
      TH2F* h_nHitsGLB_nHitsSTA;           

      // Output definition
      TFile *file_out;

};




////////
//////// -- Constructor
////////
standaloneDuplicates::standaloneDuplicates(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");

   parameters = iConfig;
   
   muonToken_  = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));
   genToken_   = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genCollection"));

   h_nSTA = new TH1F("h_nSTA", "Number of StandAlone muons; Number of events; Counts", 10, 0, 10);
   h_STA_pt = new TH1F("h_STA_pt", "StandAlone muon p_{T} (GeV); Number of muons; Counts", 100, 0, 100);
   h_STA_phi = new TH1F("h_STA_phi", "StandAlone muon #phi; Number of muons; Counts", 100, -3.14, 3.14);
   h_STA_eta = new TH1F("h_STA_pt", "StandAlone muon #eta; Number of muons; Counts", 100, -3, 3);
   h_nPairs = new TH1F("h_nPairs", "Number of STA pairs (common segment); Number of pairs; Counts", 3, 0, 3);
   h_Nextra = new TH1F("h_Nextra", "Number of StandAlone pairs (same TrackExtra); Number of pairs; Counts", 3, 0, 3);
   h_Nshared = new TH1F("h_Nshared", "Number of StandAlone pairs (common hit); Number of pairs; Counts", 3, 0, 3);
   h_nCommonSegments = new TH1F("h_nCommonSegments", "; Number of shared segments; Number of pairs", 10, 0, 10);
   h_nCommonHits = new TH1F("h_nCommonHits", "; Number of shared hits; Number of pairs", 25, 0, 25);
   h_nCommonSegments_nCommonHits = new TH2F("h_nCommonSegments_nCommonHits", "; Number of shared segments; Number of shared hits", 10, 0, 10, 25, 0, 25);
   h_nGLBsInPair = new TH1F("h_nGLBsInPair", "Number of GLB muons in pair; Counts", 3, 0, 3);
   h_ptGLB_ptSTA = new TH2F("h_ptGLB_ptSTA", ";GLB p_{T} (GeV); STA (only) p_{T} (GeV)", 40, 0, 80, 40, 0, 80);
   h_etaGLB_etaSTA = new TH2F("h_etaGLB_etaSTA2", ";GLB #eta; STA (only) #eta", 40, -2.4, 2.4, 40, -2.4, 2.4);
   h_phiGLB_phiSTA = new TH2F("h_phiGLB_phiSTA", ";GLB #phi; STA (only) #phi", 40, -3.14, 3.14, 40, -3.14, 3.14);
   h_nDTGLB_nDTSTA = new TH2F("h_nDTGLB_nDTSTA", ";GLB N_{DT}; STA (only) N_{DT}", 30, 0, 30, 30, 0, 30);
   h_nCSCGLB_nCSCSTA = new TH2F("h_nCSCGLB_nCSCSTA", ";GLB N_{CSC}; STA (only) N_{CSC}", 20, 0, 20, 20, 0, 20);
   h_nHitsGLB_nHitsSTA = new TH2F("h_nHitsGLB_nHitsSTA", ";GLB N_{Hits}; STA (only) N_{Hits}", 30, 0, 30, 30, 0, 30);
}



////////
//////// -- Destructor
////////
standaloneDuplicates::~standaloneDuplicates()
{

}

////////
//////// -- BeginJob
////////
void standaloneDuplicates::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  output_filename = parameters.getParameter<std::string>("nameOfOutput");
  file_out = new TFile(output_filename.c_str(), "RECREATE");

}

////////
//////// -- EndJob
////////
void standaloneDuplicates::endJob()
{
  std::cout << "End Job" << std::endl;

  file_out->cd();
  h_nSTA->Write();
  h_STA_pt->Write();
  h_STA_eta->Write();
  h_STA_phi->Write();
  h_nPairs->Write();
  h_nGLBsInPair->Write();
  h_nCommonSegments->Write();
  h_nCommonHits->Write();
  h_nCommonSegments_nCommonHits->Write();
  h_ptGLB_ptSTA->Write();
  h_etaGLB_etaSTA->Write();
  h_phiGLB_phiSTA->Write();
  h_nDTGLB_nDTSTA->Write();
  h_nCSCGLB_nCSCSTA->Write();
  h_nHitsGLB_nHitsSTA->Write();
  file_out->Close();

}


////////
//////// -- fillDescriptions
////////
void standaloneDuplicates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


////////
//////// -- Analyze
////////
void standaloneDuplicates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   iEvent.getByToken(muonToken_, muonCollection_);
   iEvent.getByToken(genToken_, genCollection_);

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
   

   // Generated muons
   std::vector<reco::GenParticle> genMuons;
   for (const auto& gen : *genCollection_) {
     if (fabs(gen.pdgId()) != 13)
       continue;
     if (gen.status() != 1)
       continue;
     genMuons.push_back(gen);
   }

   //std::cout << "Number of generated muons: " << genMuons.size() << std::endl; 

   // Reconstructed muons
   //std::vector<reco::Muon> muons;
   std::vector<int> iMu; // index of STA muons
   int nSTA = 0;
   for (size_t i = 0; i < muonCollection_->size(); i++) {
     const reco::Muon &muon = (*muonCollection_)[i];
     if (muon.isStandAloneMuon()) {
       iMu.push_back(i);
       nSTA++;
       h_STA_pt->Fill(muon.outerTrack()->pt());
       h_STA_eta->Fill(muon.outerTrack()->eta());
       h_STA_phi->Fill(muon.outerTrack()->phi());
     }
   }
   h_nSTA->Fill(nSTA);


   // Study duplicates
   std::vector<std::pair<reco::Muon, reco::Muon>> duplicates; // real GLB : clone STA
   std::vector<std::pair<int, int>> indexes; // real idx : clone idx
   std::vector<std::pair<int, int>> isGLB; // isGLB : isGLB
   std::vector<int> nSharedSegments;
   std::vector<int> nSharedHits;
   int Nextra = 0; // Number of pairs sharing extra track information
   int Nshared = 0; // Number of pairs sharing at least one hit

   for (unsigned int i = 0; i < iMu.size(); i++) {
     for (unsigned int j = i+1; j < iMu.size(); j++) {
       if (i==j)
         continue;

       const reco::Muon & mu_i = (*muonCollection_)[iMu.at(i)];
       const reco::Muon & mu_j = (*muonCollection_)[iMu.at(j)];

       if ( (*mu_i.outerTrack()).extra() == (*mu_j.outerTrack()).extra() )
         Nextra++;

       if (commonTrackHits( *mu_i.outerTrack(), *mu_j.outerTrack()) > 0) {

         Nshared++;

         duplicates.emplace_back( std::make_pair(mu_i, mu_j) );
         indexes.emplace_back( std::make_pair(i, j) );
         nSharedSegments.push_back(commonSegments(mu_i, mu_j));
         nSharedHits.push_back(commonTrackHits( *mu_i.outerTrack(), *mu_j.outerTrack()));
         
         h_nCommonSegments->Fill(commonSegments(mu_i, mu_j));
         h_nCommonHits->Fill(commonTrackHits( *mu_i.outerTrack(), *mu_j.outerTrack()));
         h_nCommonSegments_nCommonHits->Fill(commonSegments(mu_i, mu_j), commonTrackHits( *mu_i.outerTrack(), *mu_j.outerTrack()));

         // Evaluate if muons are Global
         reco::Muon globalMu;
         reco::Muon staMu;

         if (mu_i.isGlobalMuon() && mu_j.isGlobalMuon())
           h_nGLBsInPair->Fill(2);

         if ((mu_i.isGlobalMuon() && !mu_j.isGlobalMuon()) || (!mu_i.isGlobalMuon() && mu_j.isGlobalMuon())) {
           h_nGLBsInPair->Fill(1);
           if (mu_i.isGlobalMuon() && !mu_j.isGlobalMuon()){
             globalMu = mu_i;
             staMu = mu_j;
           } else {
             globalMu = mu_j;
             staMu = mu_i;
           }

           // Fill plot pairs:
           h_ptGLB_ptSTA->Fill(globalMu.outerTrack()->pt(), staMu.outerTrack()->pt());
           h_etaGLB_etaSTA->Fill(globalMu.outerTrack()->eta(), staMu.outerTrack()->eta());
           h_phiGLB_phiSTA->Fill(globalMu.outerTrack()->phi(), staMu.outerTrack()->phi());
           h_nDTGLB_nDTSTA->Fill(globalMu.outerTrack()->hitPattern().numberOfValidMuonDTHits(), staMu.outerTrack()->hitPattern().numberOfValidMuonDTHits());
           h_nCSCGLB_nCSCSTA->Fill(globalMu.outerTrack()->hitPattern().numberOfValidMuonCSCHits(), staMu.outerTrack()->hitPattern().numberOfValidMuonCSCHits());
           h_nHitsGLB_nHitsSTA->Fill(globalMu.outerTrack()->hitPattern().numberOfValidHits(), staMu.outerTrack()->hitPattern().numberOfValidHits());

         }

         if (!mu_i.isGlobalMuon() && !mu_j.isGlobalMuon())
           h_nGLBsInPair->Fill(0);

         /*
         int nHits = commonTrackHits( *mu_i.outerTrack(), *mu_j.outerTrack());
         std::cout << "Indexes: " << i << " " << j << std::endl;
         std::cout << "Eta: " << mu_i.eta() << " " << mu_j.eta() << std::endl;
         std::cout << "Phi: " << mu_i.phi() << " " << mu_j.phi() << std::endl;
         std::cout << "nHits: " << nHits << std::endl;
         */
       }
     }
   }

   h_nPairs->Fill(duplicates.size());
   h_Nextra->Fill(Nextra);
   h_Nshared->Fill(Nshared);
}


DEFINE_FWK_MODULE(standaloneDuplicates);
