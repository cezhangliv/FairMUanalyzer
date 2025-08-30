#include "FairMUanalyzer.h"
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <iostream>


void FairMUanalyzer::AnalyzeTRK() {

    goldenevents_ = 0;

    Long64_t N = runN_>0?runN_:cbmsim_->GetEntries();

    int TGT1 = tgt_==0?1:0;
    int TGT2 = tgt_==1?1:0;
    if(!TGT1 && !TGT2){std::cout<<"wrong set up on TGT! return."<<std::endl;return;}

    bool MF= mf_?true:false;

    std::cout << "Processing " << N << Form(" ** interaction tgt%d ** events with MF tag %s...",tgt_, mf_?"ON":"OFF") << std::endl;
    
    int igraph = 0;
    for (Long64_t i = 0; i < N; ++i) {

        if (i % (N / 10) == 0 || i == N - 1) {double progress = 100.0 * i / N;printf("Processing: %.1f%% (%lld/%lld)\n", progress, i, N);}
        case_counts["Total"]++;

        cbmsim_->GetEntry(i);
        const auto& tracks = reco_->reconstructedTracks();
        const auto& hits = reco_->reconstructedHits();
        const auto& bestvtx = reco_->bestVertex();

        int n_muons = 0;
        
        std::vector<const MUonERecoOutputTrackAnalysis*> muon_tracks;
        for (auto const& track : tracks) {
            if (track.isMuon() ) {
                n_muons++;
                muon_tracks.push_back(&track);
            }
        }
        h_isMuon->Fill(n_muons);
        h_Ntracks->Fill(tracks.size());

        int nhits_zcut = 0;
        int nhits_sec0=0; 
        int nhits_sec1=0;
        int nhits_sec2=0;

        for (auto const& hit : hits) {
            if (hit.z() > 1000) {
                nhits_zcut++;
            }

            if (hit.stationID()==0)nhits_sec0++;
            if (hit.stationID()==1)nhits_sec1++;
            if (hit.stationID()==2)nhits_sec2++;

        }
        h_hits_zcut->Fill(nhits_zcut);

        if (n_muons >= 1 && n_muons <= 4) {
            for (auto const* track : muon_tracks) {
                int nhit_zcut = 0;
                int trk_muonID = track->muonId();

                for (auto const& hit : hits) {
                    if (hit.z() > 1000) {
                        
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

        
        //golden muon step #1: N tracks

        if ( (TGT2 && (useTightTrackCutTgt2_ ? tracks.size() == 4 : tracks.size() >= 4)) 
            || 
            (TGT1 && (useTightTrackCutTgt1_?tracks.size() == 5: tracks.size() >= 3))
            ) {

            bool isGolden = true;
            std::set<int> sectors;
            std::set<int> sectors01;

            int ntrk_sec0=0; 
            int ntrk_sec1=0;
            int ntrk_sec2=0;

            for (auto const& track : tracks) {
                std::set<int> modules;
                for (auto const& h : track.hits()) {
                    modules.insert(h.moduleID());
                }
                if (modules.size() != 6   && (TGT1) && track.sector()<2 ) {
                    //golden muon step #2: 1 hit/module
                    isGolden = false;
                    break;
                }
                else if (modules.size() != 6   && (TGT2 && useTightTrackCutTgt2_)  ) {
                    //golden muon step #2: 1 hit/module
                    isGolden = false;
                    break;
                }
                else if (modules.size() < 5   && (TGT2 && !useTightTrackCutTgt2_)  ) {
                    //golden muon step #2: 1 hit/module
                    isGolden = false;
                    break;
                }

                //golden muon step #3: reduced chi2
                if(track.chi2perDegreeOfFreedom()>=2)isGolden = false;

                sectors.insert(track.sector());
                if(track.sector()<2)sectors01.insert(track.sector());

                if(track.sector()==0)ntrk_sec0++;
                if(track.sector()==1)ntrk_sec1++;
                if(track.sector()==2)ntrk_sec2++;

            }

            //golden muon step #4: 1/2 tracks/station, and nhit_giovanni cut
            //note we don't require muonid at this moment, just (1track + 2track + (1+) tracks) or (1track + 1track + 2 tracks) signature is ok
            //we have to do this way because we  need those 'wrong' events (even if we know by MF) to for the next step to see the count/distributions

            if(TGT1
                && 
                ( sectors01.size() != 2 || !(ntrk_sec0==1 && ntrk_sec1==2) || (nhits_sec1 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
                )isGolden = false;
                

            if(TGT2 
                && useTightTrackCutTgt2_
                && 
                ( sectors.size() != 3 || !(ntrk_sec0==1 && ntrk_sec1==1 && ntrk_sec2==2) || (nhits_sec2 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
                )isGolden = false;
            else if(TGT2 
                && !useTightTrackCutTgt2_
                && 
                //( sectors.size() != 3 || !(ntrk_sec0==1 && ntrk_sec1==1 && ntrk_sec2>=2) || (nhits_sec2 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
                ( sectors.size() != 3 || !(ntrk_sec0==1 && ntrk_sec1>=1 && ntrk_sec2>=2) || (nhits_sec2 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
                )isGolden = false;

            if (isGolden) {

                goldenevents_++;
                case_counts["golden"]++;

                //event selection #1: within target
                int intgt=0;
                if(TGT2 && bestvtx.zPositionFit()<=z_tgt2_+2 && bestvtx.zPositionFit()>=z_tgt2_-2)intgt=1;
                if(TGT1 && bestvtx.zPositionFit()<=z_tgt1_+2 && bestvtx.zPositionFit()>=z_tgt1_-2)intgt=1;
                //Elastic step #1: tgt position
                if(!intgt)continue;

                h_vtxchi2->Fill(bestvtx.chi2perDegreeOfFreedom());
                //Elastic step #2 (optional): bestvtx chi2perDOF
                //if(bestvtx.chi2perDegreeOfFreedom()>10)continue;

                int sec0=0; 
                int sec1=0;
                int sec2=0;

                int sec1muon=0; 
                int sec1e=0;
                int sec2muon=0;
                int sec2e=0;

                std::vector<TVector3> in; in.reserve(12);
                std::vector<TVector3> out; out.reserve(12);
                std::vector<TVector3> oute; oute.reserve(12);
                std::vector<TVector3> outmuon; outmuon.reserve(12);

                //Elastic step #3: aco (following)


                if(TGT2 && !useTightTrackCutTgt2_){

                    //if( abs(bestvtx.modifiedAcoplanarity())>0.4e-3 || bestvtx.chi2perDegreeOfFreedom()>3 )continue;//0.4 rad
                    if( abs(bestvtx.modifiedAcoplanarity())>0.4 )continue;//not very much difference between 0.4 and 0.4e-3. see my slides 250722
                    
                    case_counts["t1mem"]++;
                    case_counts["t1all"]++;
                    
                    //case_h2d["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d["t1mem"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_g2d["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_g2d["t1mem"]->SetPoint(case_g2d["t1mem"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());

                    h_2d->Fill(bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    g_2d->SetPoint(g_2d->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    
                    if(bestvtx.electronTheta()<=intersecX_){
                        case_counts["t1me<m"]++;
                        case_h2d["t1me<m"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                        case_g2d["t1me<m"]->SetPoint(case_g2d["t1me<m"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());
                    }
                    continue; // so skip the rest part using useTightTrackCutTgt2_ Ntrack==2
                }

                if(TGT2 && useTightTrackCutTgt2_){

                    if(abs(bestvtx.modifiedAcoplanarity()<0.4) )h_2d_bstvtx->Fill(bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    //if(abs(bestvtx.modifiedAcoplanarity()<0.3) )h_2d_bstvtx->Fill(bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    
                    h_2d_bstvtx->Fill(bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    g_2d_bstvtx->SetPoint(g_2d_bstvtx->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta()); 
                }
                
                for(int j=0; j<tracks.size();j++)
                {
                    if(TGT1)continue;

                    if(tracks.at(j).sector()==1) {
                        TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); in.push_back(v);

                        //Eugenia's cut https://indico.cern.ch/event/1476217/contributions/6217032/attachments/2962101/5210167/tesi_phd_weekly.pdf
                        //if(v.Theta()>4e-3)continue;
                        sec1++;
                    }
                    if(tracks.at(j).sector()==2 && MF) {
                        if(tracks.at(j).isMuon()){sec2muon++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); outmuon.push_back(v);}
                        else {sec2e++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); oute.push_back(v);}
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
                        //if(v.Theta()>4e-3 )continue;

                        sec0++;
                    }
                    if(tracks.at(j).sector()==1 && MF) {
                        if(tracks.at(j).isMuon()){sec1muon++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); outmuon.push_back(v);}
                        else {sec1e++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); oute.push_back(v);}
                    }    
                    else if(tracks.at(j).sector()==1){
                        sec1++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out.push_back(v);
                    }

                }

                double angle_e, angle_mu;
                double angle0, angle1;
                double aco = 1000;

                //case1: with MF
                if (sec1 == 1 && sec2e == 1 && sec2muon == 1 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(outmuon.at(0)); 
                    aco=acoplanarity(in.at(0),oute.at(0),outmuon.at(0)); 
                    
                    if( abs(aco)>0.4)continue;//0.4 rad
                    //if( abs(aco)>0.3)continue;//0.3 rad
                    
                    case_counts["t1all"]++;
                    case_counts["t1mem"]++;
                    case_h2d["t1all"]->Fill(angle_e,angle_mu);
                    case_h2d["t1mem"]->Fill(angle_e,angle_mu);
                    case_h2d_bstvtx["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d_bstvtx["t1mem"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_e,angle_mu);
                    case_g2d["t1mem"]->SetPoint(case_g2d["t1mem"]->GetN(),angle_e,angle_mu);
                    case_g2d_bstvtx["t1all"]->SetPoint(case_g2d_bstvtx["t1all"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_g2d_bstvtx["t1mem"]->SetPoint(case_g2d_bstvtx["t1mem"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());

                    h_2d->Fill(angle_e,angle_mu); 
                    if(angle_e<=intersecX_){
                        case_counts["t1me<m"]++;
                        case_h2d["t1me<m"]->Fill(angle_e,angle_mu);
                        case_h2d_bstvtx["t1me<m"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                        case_g2d["t1me<m"]->SetPoint(case_g2d["t1me<m"]->GetN(),angle_e,angle_mu);
                        case_g2d_bstvtx["t1me<m"]->SetPoint(case_g2d_bstvtx["t1me<m"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());
                    }
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                    //case_h2d_bstvtx["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                }    
                if (sec1 == 1 && sec2e == 2 && sec2muon == 0 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(oute.at(1));    
                    aco=acoplanarity(in.at(0),oute.at(0),oute.at(1)); 
                    
                    if( abs(aco)>0.4)continue;//0.4 rad
                    //if( abs(aco)>0.3)continue;//0.3 rad
                    
                    case_counts["t1mee"]++; 
                    case_counts["t1all"]++; 
                    if(angle_e>angle_mu){
                        case_h2d["t1mee"]->Fill(angle_e,angle_mu);case_h2d["t1all"]->Fill(angle_e,angle_mu);
                        case_g2d["t1mee"]->SetPoint(case_g2d["t1mee"]->GetN(),angle_e,angle_mu);
                        case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t1mee"]->Fill(angle_mu,angle_e);case_h2d["t1all"]->Fill(angle_mu,angle_e);
                        case_g2d["t1mee"]->SetPoint(case_g2d["t1mee"]->GetN(), angle_mu,angle_e);
                        case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(), angle_mu,angle_e);
                    }
                    
                    case_h2d_bstvtx["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d_bstvtx["t1mee"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    case_g2d_bstvtx["t1all"]->SetPoint(case_g2d_bstvtx["t1all"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_g2d_bstvtx["t1mee"]->SetPoint(case_g2d_bstvtx["t1mee"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                    //case_h2d_bstvtx["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                }    
                if (sec1 == 1 && sec2e == 0 && sec2muon == 2 && MF){ 
                    angle_e=in.at(0).Angle(outmuon.at(0)); 
                    angle_mu=in.at(0).Angle(outmuon.at(1)); 
                    aco=acoplanarity(in.at(0),outmuon.at(0),outmuon.at(1)); 
                    
                    if( abs(aco)>0.4)continue;//0.4 rad
                    //if( abs(aco)>0.3)continue;//0.3 rad
                    
                    case_counts["t1mmm"]++; 
                    case_counts["t1all"]++; 
                    
                    if(angle_e>angle_mu){
                        case_h2d["t1mmm"]->Fill(angle_e,angle_mu);case_h2d["t1all"]->Fill(angle_e,angle_mu);
                        case_g2d["t1mmm"]->SetPoint(case_g2d["t1mmm"]->GetN(),angle_e,angle_mu);
                        case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t1mmm"]->Fill(angle_mu,angle_e);case_h2d["t1all"]->Fill(angle_mu,angle_e);
                        case_g2d["t1mmm"]->SetPoint(case_g2d["t1mmm"]->GetN(),angle_mu,angle_e);
                        case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_mu,angle_e);
                    }
                    
                    case_h2d_bstvtx["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d_bstvtx["t1mmm"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    case_g2d_bstvtx["t1all"]->SetPoint(case_g2d_bstvtx["t1all"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_g2d_bstvtx["t1mmm"]->SetPoint(case_g2d_bstvtx["t1mmm"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                    //case_h2d_bstvtx["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                }    
                if (sec0 == 1 && sec1e == 1 && sec1muon == 1 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(outmuon.at(0)); 
                    aco=acoplanarity(in.at(0),oute.at(0),outmuon.at(0)); 
                    
                    if( abs(aco)>0.4)continue;//0.4 rad
                    
                    case_counts["t0mem"]++; 
                    case_counts["t0all"]++; 
                    case_h2d["t0mem"]->Fill(angle_e,angle_mu); 
                    case_h2d["t0all"]->Fill(angle_e,angle_mu); 
                    case_g2d["t0mem"]->SetPoint(case_g2d["t0mem"]->GetN(), angle_e,angle_mu); 
                    case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_e,angle_mu); 
                    
                    h_2d->Fill(angle_e,angle_mu); 
                    g_2d->SetPoint(g_2d->GetN(), angle_e,angle_mu); 
                    
                    if(angle_e<=intersecX_){
                        case_counts["t0me<m"]++;
                        case_h2d["t0me<m"]->Fill(angle_e,angle_mu); 
                        case_g2d["t0me<m"]->SetPoint(case_g2d["t0me<m"]->GetN(),angle_e,angle_mu); 
                    }
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                }    
                if (sec0 == 1 && sec1e == 2 && sec1muon == 0 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(oute.at(1));    
                    aco=acoplanarity(in.at(0),oute.at(0),oute.at(1)); 
                    
                    if( abs(aco)>0.4)continue;//0.4 rad
                    
                    case_counts["t0mee"]++; 
                    case_counts["t0all"]++; 
                    if(angle_e>angle_mu){
                        case_h2d["t0mee"]->Fill(angle_e,angle_mu);case_h2d["t0all"]->Fill(angle_e,angle_mu);
                        case_g2d["t0mee"]->SetPoint(case_g2d["t0mee"]->GetN(), angle_e,angle_mu);
                        case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t0mee"]->Fill(angle_mu,angle_e);case_h2d["t0all"]->Fill(angle_mu,angle_e);
                        case_g2d["t0mee"]->SetPoint(case_g2d["t0mee"]->GetN(), angle_mu,angle_e);
                        case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_mu,angle_e);
                    }
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                }    
                if (sec0 == 1 && sec1e == 0 && sec1muon == 2 && MF){ 
                    angle_e=in.at(0).Angle(outmuon.at(0)); 
                    angle_mu=in.at(0).Angle(outmuon.at(1)); 
                    aco=acoplanarity(in.at(0),outmuon.at(0),outmuon.at(1)); 
                    
                    if( abs(aco)>0.4)continue;//0.4 rad
                    
                    case_counts["t0mmm"]++; 
                    case_counts["t0all"]++; 
                    
                    if(angle_e>angle_mu){
                        case_h2d["t0mmm"]->Fill(angle_e,angle_mu);case_h2d["t0all"]->Fill(angle_e,angle_mu);
                        case_g2d["t0mmm"]->SetPoint(case_g2d["t0mmm"]->GetN(), angle_e,angle_mu);
                        case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t0mmm"]->Fill(angle_mu,angle_e);case_h2d["t0all"]->Fill(angle_mu,angle_e);
                        case_g2d["t0mmm"]->SetPoint(case_g2d["t0mmm"]->GetN(), angle_mu,angle_e);
                        case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_mu,angle_e);
                    }
                    
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                }    

                

                //case2: w/o MF, tgt1/2

                if( !MF && 
                    ( (sec1==1 && sec2==2) || (sec0==1 && sec1==2 ) )
                   ) 
                {

                    /*
                    double angle0=in.at(0).Angle(out.at(0)); 
                    double angle1=in.at(0).Angle(out.at(1)); 
                    //event selection #2: acoplanarity
                    double dotProduct_v = out.at(0).Dot(out.at(1));
                    TVector3 crossProduct_v = out.at(0).Cross(out.at(1));
                    double T_v = in.at(0).Dot(crossProduct_v);
                    TVector3 im_v= in.at(0).Cross(out.at(0));
                    TVector3 ie_v= in.at(0).Cross(out.at(1));
                    T_v = T_v>0? 1:-1;
                    double acoplanarity_v= T_v*(TMath::Pi() - acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));
                    */
                    
                    angle0=in.at(0).Angle(out.at(0)); 
                    angle1=in.at(0).Angle(out.at(1)); 
                    aco=acoplanarity(in.at(0),out.at(0),out.at(1));
                    if( abs(aco)>0.4)continue;//0.4 rad
            
                    //if(tracks.size()!=3 || (angle0>angle1 && angle0>0.032) || (angle0<angle1 && angle1>0.032) || (angle0>angle1 && angle1<0.0002) || (angle0<angle1 && angle0<0.0002) )continue;
                    //flag_good_event = 1;

                    if(angle0>angle1) {h_2d->Fill(angle0,angle1);g_2d->SetPoint(g_2d->GetN(),angle0,angle1);}
                    else {h_2d->Fill(angle1,angle0);g_2d->SetPoint(g_2d->GetN(),angle1,angle0);}
                    
                }

            }//isGolden

        }// simple N tracks cut
    
    }//event loop

    std::cout<<"goldenevents: "<<goldenevents_<<"/"<<N<<std::endl;
}
