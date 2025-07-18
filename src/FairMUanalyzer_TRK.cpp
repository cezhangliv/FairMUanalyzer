#include "FairMUanalyzer.h"
#include <set>
#include <algorithm>
#include <iostream>

int TGT1 = 0;
int TGT2 = 1;
bool MF=false;


void FairMUanalyzer::AnalyzeTRK() {

    goldenevents_ = 0;

    Long64_t N = cbmsim_->GetEntries();
    //Long64_t N = 100;
    std::cout << "Processing " << N << " events..." << std::endl;
    for (Long64_t i = 0; i < N; ++i) {
        cbmsim_->GetEntry(i);
        const auto& tracks = reco_->reconstructedTracks();
        const auto& hits = reco_->reconstructedHits();
        const auto& bestvtx = reco_->bestVertex();

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

        if ( (TGT2 && tracks.size() >= 4) || (TGT1 && tracks.size() >= 3) ) {

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

            //if (sectors.size() != 4) isGolden = false;

            //std::cout<<__LINE__<<std::endl;

            if(TGT2 && bestvtx.zPositionFit()<770)continue;
            if(TGT1 && bestvtx.zPositionFit()<670)continue;

            if (isGolden) {

                goldenevents_++;
                //std::cout<<__LINE__<<std::endl;

                for (auto const& track : tracks) {
                    h_goldenMuon_isMuon[track.sector()]->Fill(track.isMuon());
                }

                std::vector<TVector3> in; in.reserve(12);
                std::vector<TVector3> out; out.reserve(12);
                std::vector<TVector3> out_e; out_e.reserve(12);
                std::vector<TVector3> outmuon; outmuon.reserve(12);

                int sec0=0; 
                int sec1=0;
                int sec2=0;
                
                int sec1e=0;
                int sec1muon=0;
                int sec2e=0;
                int sec2muon=0;


                for(int j=0; j<tracks.size();j++)
                {
                    if(TGT1)continue;

                    if(tracks.at(j).sector()==1) {
                        TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); in.push_back(v);
                        
                        //Eugenia's cut https://indico.cern.ch/event/1476217/contributions/6217032/attachments/2962101/5210167/tesi_phd_weekly.pdf
                        if(v.Theta()>4e-3 || tracks.at(j).chi2perDegreeOfFreedom()>=2 )continue;

                        sec1++;
                    }
                    if(tracks.at(j).sector()==2 && MF) {
                        if(tracks.at(j).isMuon()){sec2muon++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); outmuon.push_back(v);}
                        else {sec2e++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out_e.push_back(v);}
                    }    
                    else if(tracks.at(j).sector()==2){
                        sec2++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out.push_back(v);
                        

                    }

                }

                for(int j=0; j<tracks.size();j++)
                {
                    if(TGT2)continue;
                    
                    if(tracks.at(j).sector()==0) {
                        TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); in.push_back(v);
                        
                        //Eugenia's cut https://indico.cern.ch/event/1476217/contributions/6217032/attachments/2962101/5210167/tesi_phd_weekly.pdf
                        if(v.Theta()>4e-3 || tracks.at(j).chi2perDegreeOfFreedom()>=2 )continue;

                        sec0++;
                    }
                    if(tracks.at(j).sector()==1 && MF) {
                        if(tracks.at(j).isMuon()){sec1muon++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); outmuon.push_back(v);}
                        else {sec1e++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out_e.push_back(v);}
                    }    
                    else if(tracks.at(j).sector()==1){
                        sec1++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out.push_back(v);
                        

                    }

                }
                //case1: with MF tgt1/2
                if( MF && 
                    ( (sec1==1 && sec2e==1 && sec2muon==1) || (sec0==1 && sec1e==1 && sec1muon==1) )
                   ) 
                {

                    //std::cout<<__LINE__<<std::endl;
                    //double angle0=in.at(0).Angle(out.at(0)); 
                    //double angle1=in.at(0).Angle(out.at(1)); 
                    
                    double angle_e=in.at(0).Angle(out_e.at(0)); 
                    double angle_mu=in.at(0).Angle(outmuon.at(0)); 


                    double dotProduct_v = outmuon.at(0).Dot(out_e.at(0));
                    TVector3 crossProduct_v = outmuon.at(0).Cross(out_e.at(0));

                    double T_v = in.at(0).Dot(crossProduct_v);
                    TVector3 im_v= in.at(0).Cross(outmuon.at(0));
                    TVector3 ie_v= in.at(0).Cross(out_e.at(0));
                    T_v = T_v>0? 1:-1;
                    double acoplanarity_v= T_v*(TMath::Pi() - acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

                    if( abs(acoplanarity_v)>0.4)continue;//0.4 rad
            
                    //if(tracks.size()!=3 || (angle0>angle1 && angle0>0.032) || (angle0<angle1 && angle1>0.032) || (angle0>angle1 && angle1<0.0002) || (angle0<angle1 && angle0<0.0002) )continue;
                    //flag_good_event = 1;

                    //if(angle0>angle1) h_2d->Fill(angle0,angle1);
                    //else h_2d->Fill(angle1,angle0);
                    h_2d->Fill(angle_e,angle_mu);
                }

                //case2: w/o MF, tgt1/2
                if( !MF && 
                    ( (sec1==1 && sec2==2) || (sec0==1 && sec1==2 ) )
                   ) 
                {

                    //std::cout<<__LINE__<<std::endl;
                    double angle0=in.at(0).Angle(out.at(0)); 
                    double angle1=in.at(0).Angle(out.at(1)); 



                    double dotProduct_v = out.at(0).Dot(out.at(1));
                    TVector3 crossProduct_v = out.at(0).Cross(out.at(1));

                    double T_v = in.at(0).Dot(crossProduct_v);
                    TVector3 im_v= in.at(0).Cross(out.at(0));
                    TVector3 ie_v= in.at(0).Cross(out.at(1));
                    T_v = T_v>0? 1:-1;
                    double acoplanarity_v= T_v*(TMath::Pi() - acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

                    if( abs(acoplanarity_v)>0.4)continue;//0.4 rad
            
                    //if(tracks.size()!=3 || (angle0>angle1 && angle0>0.032) || (angle0<angle1 && angle1>0.032) || (angle0>angle1 && angle1<0.0002) || (angle0<angle1 && angle0<0.0002) )continue;
                    //flag_good_event = 1;

                    if(angle0>angle1) h_2d->Fill(angle0,angle1);
                    else h_2d->Fill(angle1,angle0);
                    
                }





                /*
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
                */
            }
        }
    }

    std::cout<<"goldenevents: "<<goldenevents_<<std::endl;
}