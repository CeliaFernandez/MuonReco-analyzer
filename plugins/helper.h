//
// Original Author:
//                Celia Fernandez Madrazo
//
// helper functions
//

#ifndef MuonReco_analysis_plugins_helper
#define MuonReco_analysis_plugins_helper

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"


// Common segments between two muons accessed with embedded matches

template <typename MUO>
inline int commonSegments(const MUO muon1, const MUO muon2)
{
  int nmatches = 0;
  for (auto & chamber1 : muon1.matches()) {
    for (auto & chamber2 : muon2.matches()) {
      if (chamber1.id.rawId() != chamber2.id.rawId())
        continue;
      if (chamber1.id.subdetId() != chamber2.id.subdetId())
        continue;
      if (chamber1.station() != chamber2.station())
        continue;
      for (auto & segment1 : chamber1.segmentMatches) {
        for (auto & segment2 : chamber2.segmentMatches) {
          if (fabs(segment1.x - segment2.x) < 1e-6 && fabs(segment1.y - segment2.y) < 1e-6) {
            nmatches++;
            break;
          }
        }
      }
    }
  }

  return nmatches;

}


// Common hits between two muons accessed with embedded matches

inline int commonTrackHits(const reco::Track track1, const reco::Track track2)
{
  int nsharedhits = 0;
  for (auto& hit1 : track1.recHits()){
    if (!hit1->isValid())
      continue;
    DetId id1 = hit1->geographicalId();
    //if (id1.det() != DetId::Muon)
    //  continue;
    for (auto& hit2 : track2.recHits()){
      if (!hit2->isValid())
        continue;
      DetId id2 = hit2->geographicalId();
      if (id1.rawId() != id2.rawId())
        continue;
      if (fabs(hit1->localPosition().x() - hit2->localPosition().x()) < 1e-6 && fabs(hit1->localPosition().y() - hit2->localPosition().y()) < 1e-6) {
        nsharedhits++;
        break;
      }
    }
  }

  return nsharedhits;

}

#endif


