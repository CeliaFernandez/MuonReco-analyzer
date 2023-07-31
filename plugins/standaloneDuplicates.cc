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

      // Track collections
      edm::EDGetTokenT<edm::View<reco::Track> > staToken_;     // standAloneMuons
      edm::Handle<edm::View<reco::Track> > staTracks_;
      edm::EDGetTokenT<edm::View<reco::Track> > staVtxToken_;  // standAloneMuons::UpdatedAtVtx
      edm::Handle<edm::View<reco::Track> > staVtxTracks_;
      edm::EDGetTokenT<edm::View<reco::Track> > glbToken_;     // globalMuons
      edm::Handle<edm::View<reco::Track> > glbTracks_;

      // Muon collection
      edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
      edm::Handle<edm::View<reco::Muon> > muonCollection_;

      //
      // --- Variables used
      //

      // Event info
      Int_t eventId = 0;
      Int_t luminosityBlock = 0;
      Int_t run = 0;



      // Histograms definition

      // Histograms used to study the beamspot constraint, uses both reco::Track and reco::Muon objects
      TH1F* h_success;
      TH1F* h_gsta_pt;
      TH1F* h_bsta_pt;
      TH1F* h_gsta_eta;
      TH1F* h_bsta_eta;
      TH1F* h_insta_pt;
      TH1F* h_outsta_pt;
      TH1F* h_insta_eta;
      TH1F* h_outsta_eta;
      TH2F* h_insta_pt_eta;
      TH2F* h_outsta_pt_eta;
      TH2F* h_nonupd_upd_pt;
      TH2F* h_nonupd_upd_eta;
      TH1F* h_allmuons_pt;
      TH1F* h_allmuons_eta;
      TH2F* h_allmuons_pt_eta;
      TH1F* h_inmuons_pt;
      TH1F* h_inmuons_eta;
      TH2F* h_inmuons_pt_eta;
      TH1F* h_outmuons_pt;
      TH1F* h_outmuons_eta;
      TH2F* h_outmuons_pt_eta;


      // reco::Muon Histograms (used to evaluate the presence of duplicates in the reco::Muon colletion)
      TH1F* h_nSTA;
      TH1F* h_STA_pt;
      TH1F* h_STA_eta;
      TH1F* h_STA_phi;
      TH1F* h_nGLB;
      TH1F* h_GLB_pt;
      TH1F* h_GLB_eta;
      TH1F* h_GLB_phi;
      TH1F* h_nPairs;
      TH1F* h_Nextra;
      TH1F* h_Nshared;
      TH1F* h_nGLBsInPair;
      TH1F* h_nCommonSegments;
      TH2F* h_nCommonSegments_nCommonHits;
      TH1F* h_nCommonHits;
      TH2F* h_ptGLB_ptSTA;           
      TH2F* h_etaGLB_etaSTA;           
      TH2F* h_etaGLB_etaSTA_loweta;           
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
   
   muonToken_     = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));
   staToken_      = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("standAloneMuons"));
   staVtxToken_   = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("standAloneMuonsVtx"));
   glbToken_      = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("globalMuons"));

   h_success                      = new TH1F("h_success", "standAloneMuons; Type; Tracks", 4, 0, 4);
   h_gsta_pt                      = new TH1F("h_gSTA_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_bsta_pt                      = new TH1F("h_bSTA_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_gsta_eta                     = new TH1F("h_gSTA_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_bsta_eta                     = new TH1F("h_bSTA_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_insta_pt                     = new TH1F("h_inSTA_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_outsta_pt                    = new TH1F("h_outSTA_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_insta_eta                    = new TH1F("h_inSTA_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_outsta_eta                   = new TH1F("h_outSTA_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_insta_pt_eta                 = new TH2F("h_inSTA_pt_eta", "; Track p_{T} (non-updated); Track #eta (non-updated)", 100, 0, 100, 100, -3, 3);
   h_outsta_pt_eta                = new TH2F("h_outSTA_pt_eta", "; Track p_{T} (non-updated); Track #eta (non-updated)", 100, 0, 100, 100, -3, 3);
   h_nonupd_upd_pt                = new TH2F("h_upd_nonupd_pt", ";standAlone track p_{T}; standAlone track (updated) p_{T}", 50, 0, 100, 50, 0, 100);
   h_nonupd_upd_eta               = new TH2F("h_upd_nonupd_eta", ";standAlone track #eta; standAlone track (updated) #eta", 100, -3, 3, 100, -3, 3);
   h_allmuons_pt                  = new TH1F("h_allmuons_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_allmuons_eta                 = new TH1F("h_allmuons_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_allmuons_pt_eta              = new TH2F("h_allmuons_pt_eta", "; Track p_{T} (non-updated); Tracks", 100, 0, 100, 100, -3, 3);
   h_inmuons_pt                   = new TH1F("h_inmuons_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_inmuons_eta                  = new TH1F("h_inmuons_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_inmuons_pt_eta               = new TH2F("h_inmuons_pt_eta", "; Track p_{T} (non-updated); Tracks", 100, 0, 100, 100, -3, 3);
   h_outmuons_pt                  = new TH1F("h_outmuons_pt", "; Track p_{T} (non-updated); Tracks", 100, 0, 100);
   h_outmuons_eta                 = new TH1F("h_outmuons_eta", "; Track #eta (non-updated); Tracks", 100, -3, 3);
   h_outmuons_pt_eta              = new TH2F("h_outmuons_pt_eta", "; Track p_{T} (non-updated); Tracks", 100, 0, 100, 100, -3, 3);

   h_nSTA                        = new TH1F("h_nSTA", "Number of StandAlone muons; Number of events; Counts", 10, 0, 10);
   h_STA_pt                      = new TH1F("h_STA_pt", "StandAlone muon p_{T} (GeV); Number of muons; Counts", 100, 0, 100);
   h_STA_phi                     = new TH1F("h_STA_phi", "StandAlone muon #phi; Number of muons; Counts", 100, -3.14, 3.14);
   h_STA_eta                     = new TH1F("h_STA_eta", "StandAlone muon #eta; Number of muons; Counts", 100, -3, 3);
   h_nGLB                        = new TH1F("h_nGLB", "Number of Global muons; Number of events; Counts", 10, 0, 10);
   h_GLB_pt                      = new TH1F("h_GLB_pt", "Global muon p_{T} (GeV); Number of muons; Counts", 100, 0, 100);
   h_GLB_phi                     = new TH1F("h_GLB_phi", "Global muon #phi; Number of muons; Counts", 100, -3.14, 3.14);
   h_GLB_eta                     = new TH1F("h_GLB_eta", "Global muon #eta; Number of muons; Counts", 100, -3, 3);
   h_nPairs                      = new TH1F("h_nPairs", "Number of STA pairs (common segment); Number of pairs; Counts", 3, 0, 3);
   h_Nextra                      = new TH1F("h_Nextra", "Number of StandAlone pairs (same TrackExtra); Number of pairs; Counts", 3, 0, 3);
   h_Nshared                     = new TH1F("h_Nshared", "Number of StandAlone pairs (common hit); Number of pairs; Counts", 3, 0, 3);
   h_nCommonSegments             = new TH1F("h_nCommonSegments", "; Number of shared segments; Number of pairs", 10, 0, 10);
   h_nCommonHits                 = new TH1F("h_nCommonHits", "; Number of shared hits; Number of pairs", 25, 0, 25);
   h_nCommonSegments_nCommonHits = new TH2F("h_nCommonSegments_nCommonHits", "; Number of shared segments; Number of shared hits", 10, 0, 10, 25, 0, 25);
   h_nGLBsInPair                 = new TH1F("h_nGLBsInPair", "Number of GLB muons in pair; Counts", 3, 0, 3);
   h_ptGLB_ptSTA                 = new TH2F("h_ptGLB_ptSTA", ";STA (from GLB) p_{T} (GeV); STA (only) p_{T} (GeV)", 40, 0, 80, 40, 0, 80);
   h_etaGLB_etaSTA               = new TH2F("h_etaGLB_etaSTA", ";STA (from GLB) #eta; STA (only) #eta", 40, -2.4, 2.4, 40, -2.4, 2.4);
   h_etaGLB_etaSTA_loweta        = new TH2F("h_etaGLB_etaSTA_loweta", ";STA (from GLB) #eta; STA (only) #eta", 40, -0.01, 0.01, 40, -0.01, 0.01);
   h_phiGLB_phiSTA               = new TH2F("h_phiGLB_phiSTA", ";STA (from GLB) #phi; STA (only) #phi", 40, -3.14, 3.14, 40, -3.14, 3.14);
   h_nDTGLB_nDTSTA               = new TH2F("h_nDTGLB_nDTSTA", ";STA (from GLB) N_{DT}; STA (only) N_{DT}", 30, 0, 30, 30, 0, 30);
   h_nCSCGLB_nCSCSTA             = new TH2F("h_nCSCGLB_nCSCSTA", ";STA (from GLB) N_{CSC}; STA (only) N_{CSC}", 20, 0, 20, 20, 0, 20);
   h_nHitsGLB_nHitsSTA           = new TH2F("h_nHitsGLB_nHitsSTA", ";STA (from GLB) N_{Hits}; STA (only) N_{Hits}", 30, 0, 30, 30, 0, 30);
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

  h_success->Write();
  h_gsta_pt->Write();
  h_bsta_pt->Write();
  h_gsta_eta->Write();
  h_bsta_eta->Write();
  h_insta_pt->Write();
  h_outsta_pt->Write();
  h_insta_eta->Write();
  h_outsta_eta->Write();
  h_outsta_pt_eta->Write();
  h_insta_pt_eta->Write();
  h_nonupd_upd_pt->Write();
  h_nonupd_upd_eta->Write();
  h_allmuons_pt->Write();
  h_allmuons_eta->Write();
  h_allmuons_pt_eta->Write();
  h_inmuons_pt->Write();
  h_inmuons_eta->Write();
  h_inmuons_pt_eta->Write();
  h_outmuons_pt->Write();
  h_outmuons_eta->Write();
  h_outmuons_pt_eta->Write();

  h_nSTA->Write();
  h_STA_pt->Write();
  h_STA_eta->Write();
  h_STA_phi->Write();
  h_nGLB->Write();
  h_GLB_pt->Write();
  h_GLB_eta->Write();
  h_GLB_phi->Write();
  h_nPairs->Write();
  h_Nextra->Write();
  h_Nshared->Write();
  h_nGLBsInPair->Write();
  h_nCommonSegments->Write();
  h_nCommonHits->Write();
  h_nCommonSegments_nCommonHits->Write();
  h_ptGLB_ptSTA->Write();
  h_etaGLB_etaSTA->Write();
  h_etaGLB_etaSTA_loweta->Write();
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
   iEvent.getByToken(staToken_, staTracks_);
   iEvent.getByToken(staVtxToken_, staVtxTracks_);
   iEvent.getByToken(glbToken_, glbTracks_);

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

   ////   Reconstructed standaAlone tracks
   std::vector<std::pair<reco::Track, reco::Track>> pair_refits; // successfully refitted standalone and refit
   std::vector<reco::Track> failed_refits;                       // unsuccessfully refitted standalones
   std::vector<reco::Track> lost_sta;                            // unsuccessfully refitted standalones that do not enter the collection

   // Identify which standalone tracks have a good refit and which don't
   bool goodRefit;
   goodRefit = false;
   bool lostMuon;
   lostMuon = false;
   bool isGlobalMuon;
   isGlobalMuon = false;
   bool inMuons;
   inMuons = false;

   for (size_t i = 0; i < staTracks_->size(); i++) { 

     goodRefit = false;
     const reco::Track &sta = (*staTracks_)[i];
     h_success->Fill(0);

     for (size_t j = 0; j < staVtxTracks_->size(); j++) {

       const reco::Track &staVtx = (*staVtxTracks_)[j];

       if (sta.extra().get() == staVtx.extra().get()) {
         pair_refits.emplace_back( std::make_pair(sta, staVtx) );
         goodRefit = true;
         h_nonupd_upd_pt->Fill(sta.pt(), staVtx.pt());
         h_nonupd_upd_eta->Fill(sta.eta(), staVtx.eta());
         break;
       }

     }

     // Check if the sta is in muons
     inMuons = false;
     for (size_t j = 0; j < muonCollection_->size(); j++) {
       const reco::Muon &muon = (*muonCollection_)[j];
       if (!muon.isStandAloneMuon())
         continue;
       if (muon.standAloneMuon()->extra().get() == sta.extra().get()) {
         inMuons = true;
         break; 
       }
     }
     h_allmuons_pt->Fill(sta.pt());
     h_allmuons_eta->Fill(sta.eta());
     h_allmuons_pt_eta->Fill(sta.pt(), sta.eta());
     if (inMuons) {
       h_inmuons_pt->Fill(sta.pt());
       h_inmuons_eta->Fill(sta.eta());
       h_inmuons_pt_eta->Fill(sta.pt(), sta.eta());
     } else {
       h_outmuons_pt->Fill(sta.pt());
       h_outmuons_eta->Fill(sta.eta());
       h_outmuons_pt_eta->Fill(sta.pt(), sta.eta());
     }

     if (!goodRefit) {

       h_success->Fill(1);
       failed_refits.push_back(sta);
       lostMuon = true; 
       isGlobalMuon = false; 

       // Check how many of them make it into global muons
       for (size_t j = 0; j < glbTracks_->size(); j++) {
         const reco::Track &glb = (*glbTracks_)[j];
         // if (deltaR(glb, sta) < 0.15) {
         if (glb.outerDetId() == sta.outerDetId()) {
           isGlobalMuon = true;
           //std::cout << glb.outerDetId() << "\t" << sta.outerDetId() << std::endl;
           break;
         }
       }

       if (!isGlobalMuon)
         h_success->Fill(2);

       // Check how many of them make into the muons collection
       for (size_t j = 0; j < muonCollection_->size(); j++) {
         const reco::Muon &muon = (*muonCollection_)[j];
         if (!muon.isStandAloneMuon())
           continue;
         //if (muon.standAloneMuon().get() == &sta) {
         //  lostMuon = false;
         //  break; 
         //}
         if (muon.standAloneMuon()->extra().get() == sta.extra().get()) {
           lostMuon = false;
           break; 
         }
       }

       if (lostMuon) {
         h_success->Fill(3);
         lost_sta.push_back(sta);
       }

       h_bsta_pt->Fill(sta.pt());
       h_bsta_eta->Fill(sta.eta());

     } else {
       h_gsta_pt->Fill(sta.pt());
       h_gsta_eta->Fill(sta.eta());

     }

   }

   for (size_t i = 0; i < lost_sta.size(); i++) {
     const reco::Track &lost = lost_sta.at(i);
     h_outsta_pt->Fill(lost.pt());
     h_outsta_eta->Fill(lost.eta());
     h_outsta_pt_eta->Fill(lost.pt(), lost.eta());
   }

   //std::cout << "Good refits: " << pair_refits.size() << "; Bad refits: " << failed_refits.size() << "; lost sta: " << lost_sta.size() << std::endl;


   ////    Reconstructed muons
   //std::vector<reco::Muon> muons;
   std::vector<int> iMu; // index of STA muons
   int nSTA = 0;
   int nGLB = 0;
   for (size_t i = 0; i < muonCollection_->size(); i++) {
     const reco::Muon &muon = (*muonCollection_)[i];
     if (muon.isStandAloneMuon()) {
       iMu.push_back(i);
       nSTA++;
       h_STA_pt->Fill(muon.outerTrack()->pt());
       h_STA_eta->Fill(muon.outerTrack()->eta());
       h_STA_phi->Fill(muon.outerTrack()->phi());
     }
     if (muon.isGlobalMuon()) {
       nGLB++;
       h_GLB_pt->Fill(muon.globalTrack()->pt());
       h_GLB_eta->Fill(muon.globalTrack()->eta());
       h_GLB_phi->Fill(muon.globalTrack()->phi());
     }
   }
   h_nSTA->Fill(nSTA);
   h_nGLB->Fill(nGLB);


   // Study duplicates in reco::Muon collection
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
           h_etaGLB_etaSTA_loweta->Fill(globalMu.outerTrack()->eta(), staMu.outerTrack()->eta());
           h_phiGLB_phiSTA->Fill(globalMu.outerTrack()->phi(), staMu.outerTrack()->phi());
           h_nDTGLB_nDTSTA->Fill(globalMu.outerTrack()->hitPattern().numberOfValidMuonDTHits(), staMu.outerTrack()->hitPattern().numberOfValidMuonDTHits());
           h_nCSCGLB_nCSCSTA->Fill(globalMu.outerTrack()->hitPattern().numberOfValidMuonCSCHits(), staMu.outerTrack()->hitPattern().numberOfValidMuonCSCHits());
           h_nHitsGLB_nHitsSTA->Fill(globalMu.outerTrack()->hitPattern().numberOfValidHits(), staMu.outerTrack()->hitPattern().numberOfValidHits());

         }

         if (!mu_i.isGlobalMuon() && !mu_j.isGlobalMuon())
           h_nGLBsInPair->Fill(0);

       }
     }
   }

   h_nPairs->Fill(duplicates.size());
   h_Nextra->Fill(Nextra);
   h_Nshared->Fill(Nshared);
}


DEFINE_FWK_MODULE(standaloneDuplicates);
