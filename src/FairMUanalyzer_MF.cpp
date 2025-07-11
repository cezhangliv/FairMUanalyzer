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

        int n_muons = 0;
        std::vector<const MUonERecoOutputTrackAnalysis*> muon_tracks;

        for (auto const& track : tracks) {
            if (track.isMuon() && track.sector()==2 ) {
                n_muons++;
                muon_tracks.push_back(&track);
            }
        }
        h_isMuon->Fill(n_muons);
        h_Ntracks->Fill(tracks.size());

        int nhits_zcut = 0;
        for (auto const& hit : hits) {
            if (hit.z() > 1000) {
                nhits_zcut++;
            }
        }
        h_hits_zcut->Fill(nhits_zcut);

        if (n_muons >= 1 && n_muons <= 4) {
            for (auto const* track : muon_tracks) {
                int nhit_zcut = 0;
                int trk_muonID = track->muonId();

                for (auto const& hit : hits) {
                    if (hit.z() > 1000) {
                        
                        //cout<<"module ID: "<<hit.moduleID()<<endl;
                        h_hitsModuleID_zcut[n_muons]->Fill(hit.moduleID());

                        auto const& muIDs = hit.muonIds();  // vector<Int_t>
                        if (std::find(muIDs.begin(), muIDs.end(), trk_muonID) != muIDs.end()) {
                            nhit_zcut++;
                        }
                    }
                }

                h_hitsPerMuonTrack_zcut[n_muons]->Fill(nhit_zcut);
            }
        }

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

                for (auto const& track : tracks) {
                    h_goldenMuon_isMuon[track.sector()]->Fill(track.isMuon());
                }

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