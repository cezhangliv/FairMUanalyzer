#include "FairMUanalyzer.h"
#include <set>
#include <algorithm>
#include <iostream>

void FairMUanalyzer::AnalyzeMF() {
    Long64_t N = cbmsim_->GetEntries();
    std::cout << "Processing " << N << " events..." << std::endl;
    for (Long64_t i = 0; i < N; ++i) {
        cbmsim_->GetEntry(i);
        const auto& tracks = reco_->reconstructedTracks();
        const auto& hits = reco_->reconstructedHits();

        if (tracks.size() == 3 &&
            tracks[0].hits().size() == 6 &&
            tracks[1].hits().size() == 6 &&
            tracks[2].hits().size() == 6) {

            bool isGolden = true;
            std::set<int> sectors;
            for (auto const& track : tracks) {
                std::set<int> modules;
                for (auto const& h : track.hits()) {
                    modules.insert(h.moduleID());
                }
                if (modules.size() != 6) {
                    isGolden = false;
                    break;
                }
                sectors.insert(track.sector());
            }

            if (sectors.size() != 3) isGolden = false;

            if (isGolden) {
                for (int t = 0; t < 3; ++t) {
                    const auto& muonTrack = tracks[t];
                    int muID = muonTrack.muonId();
                    TVector3 p(muonTrack.xSlope(), muonTrack.ySlope(), 1.0);
                    p = p.Unit();
                    TVector3 x0(muonTrack.x0(), muonTrack.y0(), muonTrack.z0());

                    for (const auto& hit : hits) {
                        if (hit.z() < 1000) continue;
                        if (hit.moduleID() > 3 || hit.stationID() != 3) continue;
                        TVector3 pos = getXYfromHitMF(hit);
                        double dist = computeSigned2DResidualMF(p, x0, pos, hit.moduleID());

                        const auto& muIDs = hit.muonIds();
                        if (std::find(muIDs.begin(), muIDs.end(), muID) != muIDs.end()) {
                            h_residual_hitOnTrack[t]->Fill(dist);
                            h_residual_hitOnTrackModule[t][hit.moduleID()]->Fill(dist);
                            h_residual_hitAllTrackModule[t][hit.moduleID()]->Fill(dist);
                        } else {
                            h_residual_hitOffTrack[t]->Fill(dist);
                            h_residual_hitOffTrackModule[t][hit.moduleID()]->Fill(dist);
                            h_residual_hitAllTrackModule[t][hit.moduleID()]->Fill(dist);
                        }
                    }
                }
            }
        }
    }
}