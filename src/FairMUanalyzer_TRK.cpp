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
        

        cbmsim_->GetEntry(i);

        const auto& tracks = reco_->reconstructedTracks();
        const auto& hits = reco_->reconstructedHits();
        const auto& bestvtx = reco_->bestVertex();

        /// total
        case_counts["Total"]++;
        case_h1d_vertex["Total"]->Fill(bestvtx.zPositionFit());
        case_h1d_Vtxchi2["Total"]->Fill(bestvtx.chi2perDegreeOfFreedom());

        h_totalE->Fill(totalE_);
        h_clusterE->Fill(clusterE_);
        h_Ntracks->Fill(tracks.size());

        // apply cut: ECAL cluster energy <= 2/3 GeV
        //if (clusterE_ < 2.0) continue;//603
        //if (clusterE_ < 3.0) continue;//604
        //if (totalE_ < 2.0) continue;//605
        //if (totalE_ < 3.0) continue;//606

        std::vector<const MUonERecoOutputHitAnalysis*> LeftOverHits0;
        std::vector<const MUonERecoOutputHitAnalysis*> LeftOverHits1;
        std::vector<const MUonERecoOutputHitAnalysis*> LeftOverHits2;

        int nhits_MF = 0;//MF
        int nhits_sec0=0; 
        int nhits_sec1=0;
        int nhits_sec2=0;

        // MF checks
        for (auto const& hit : hits) {
            if (hit.z() > 1000) {
                nhits_MF++;
            }

            if (hit.stationID()==0)nhits_sec0++;
            if (hit.stationID()==1)nhits_sec1++;
            if (hit.stationID()==2)nhits_sec2++;

            bool found = false;

            for (const auto& track : tracks) {
                for (const auto& trkHit : track.hits()) {
                    if (trkHit.index() == hit.index()) {
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!found && hit.z()<1000) { // some hits not found because from MF
                switch (hit.stationID()) {
                    case 0:
                        LeftOverHits0.push_back(&hit);
                        break;
                    case 1:
                        LeftOverHits1.push_back(&hit);
                        break;
                    case 2:
                        LeftOverHits2.push_back(&hit);
                        break;
                    default:
                        std::cerr << "Warning: unexpected stationID " << hit.stationID()
                                  << " for hit index " << hit.index() << std::endl;
                        break;
                }
            }

        }
        h_hits_zcut->Fill(nhits_MF);

        
        //golden muon step #1: N tracks

        if(useTightTrackCutTgt2_ && tracks.size()==4)std::cout<<"good"<<endl;

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

                if (track.hits().size() != 6) {
                    isGolden = false;
                    break;
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
                //if(track.chi2perDegreeOfFreedom()>=2)isGolden = false;

                sectors.insert(track.sector());
                if(track.sector()<2)sectors01.insert(track.sector());

                if(track.sector()==0)ntrk_sec0++;
                if(track.sector()==1)ntrk_sec1++;
                if(track.sector()==2)ntrk_sec2++;

            }

            if(isGolden==true){

                case_h1d_vertex["TotalHitCut"]->Fill(bestvtx.zPositionFit());
                case_h1d_Vtxchi2["TotalHitCut"]->Fill(bestvtx.chi2perDegreeOfFreedom());
            }

            //golden muon step #4: 1/2 tracks/station, and nhit_giovanni cut
            //note we don't require muonid at this moment, just (1track + 2track + (1+) tracks) or (1track + 1track + 2 tracks) signature is ok
            //we have to do this way because we  need those 'wrong' events (even if we know by MF) to for the next step to see the count/distributions

            //h_Nhits0
            //h_Nhits1
            //h_Nhits2

            if(TGT1
                && 
                ( sectors01.size() != 2 || !(ntrk_sec0==1 && ntrk_sec1==2) ||(nhits_sec1 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
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

            if(isGolden==true){

                h_Nhits0->Fill(nhits_sec0);
                h_Nhits1->Fill(nhits_sec1);
                h_Nhits2->Fill(nhits_sec2);

                case_h1d_vertex["TotalNTrackCut"]->Fill(bestvtx.zPositionFit());
                case_h1d_Vtxchi2["TotalNTrackCut"]->Fill(bestvtx.chi2perDegreeOfFreedom());

            }

            bool LeftOverHit = false;
            
            if(TGT1 && ( ! (LeftOverHits0.size()==0 && LeftOverHits1.size()==0)  ) ){
                LeftOverHit = true;//isGolden = false;    
            }
            
            if(TGT2 
                && useTightTrackCutTgt2_
                && 
                (  !  ( LeftOverHits1.size()==0 ) ) 
                //(  !  (LeftOverHits0.size()==0 && LeftOverHits1.size()==0 && LeftOverHits2.size()==0) ) 
                ){
                LeftOverHit = true;//isGolden = false;
            }
            else if(TGT2 
                && !useTightTrackCutTgt2_
                && 
                ( !(LeftOverHits0.size()==0 ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
                )LeftOverHit = true;//isGolden = false;
            

            if(LeftOverHit){

                case_h1d_vertex["AllLeftOverhit"]->Fill(bestvtx.zPositionFit());
                case_h1d_Vtxchi2["AllLeftOverhit"]->Fill(bestvtx.chi2perDegreeOfFreedom());

                case_h1d_LeftOverhits0["AllLeftOverhit"]->Fill(LeftOverHits0.size());
                case_h1d_LeftOverhits1["AllLeftOverhit"]->Fill(LeftOverHits1.size());
                case_h1d_LeftOverhits2["AllLeftOverhit"]->Fill(LeftOverHits2.size());
                
                //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["AllLeftOverhit"]->Fill(LeftOverHits0.size());
                //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["AllLeftOverhit"]->Fill(LeftOverHits1.size());
                //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["AllLeftOverhit"]->Fill(LeftOverHits2.size());
                processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "AllLeftOverhit");
                processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "AllLeftOverhit");
                processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "AllLeftOverhit");
            }

            //if(LeftOverHit)continue;

            if (isGolden) {

                goldenevents_++;
                case_counts["golden"]++;

                case_h1d_vertex["golden"]->Fill(bestvtx.zPositionFit());
                case_h1d_Vtxchi2["golden"]->Fill(bestvtx.chi2perDegreeOfFreedom());

                case_h1d_LeftOverhits0["golden"]->Fill(LeftOverHits0.size());
                case_h1d_LeftOverhits1["golden"]->Fill(LeftOverHits1.size());
                case_h1d_LeftOverhits2["golden"]->Fill(LeftOverHits2.size());
                
                //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["golden"]->Fill(LeftOverHits0.size());
                //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["golden"]->Fill(LeftOverHits1.size());
                //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["golden"]->Fill(LeftOverHits2.size());
                processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "golden");
                processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "golden");
                processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "golden");

                //event selection #1: within target
                int intgt=0;
                if(TGT2 && bestvtx.zPositionFit()<=z_tgt2_+2 && bestvtx.zPositionFit()>=z_tgt2_-2)intgt=1;
                if(TGT1 && bestvtx.zPositionFit()<=z_tgt1_+2 && bestvtx.zPositionFit()>=z_tgt1_-2)intgt=1;
                //Elastic step #1: tgt position
                //if(!intgt)continue;

                //h_vtxchi2->Fill(bestvtx.chi2perDegreeOfFreedom());
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

                std::vector<const MUonERecoOutputTrackAnalysis*> muone_in; muone_in.reserve(12);
                std::vector<const MUonERecoOutputTrackAnalysis*> muone_out; muone_out.reserve(12);
                std::vector<const MUonERecoOutputTrackAnalysis*> muone_oute; muone_oute.reserve(12);
                std::vector<const MUonERecoOutputTrackAnalysis*> muone_outmuon; muone_outmuon.reserve(12);

                //Elastic step #3: aco (following)
                int acocut = 0;


                if(TGT2 && !useTightTrackCutTgt2_){

                    //if( abs(bestvtx.modifiedAcoplanarity())>0.4e-3 || bestvtx.chi2perDegreeOfFreedom()>3 )continue;//0.4 rad
                    if( abs(bestvtx.modifiedAcoplanarity())>0.4 && acocut)continue;//not very much difference between 0.4 and 0.4e-3. see my slides 250722
                    
                    case_counts["t1mem"]++;
                    case_counts["t1all"]++;
                    
                    //case_h2d["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d["t1mem"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_g2d["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    //case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());
                    //case_g2d["t1mem"]->SetPoint(case_g2d["t1mem"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());

                    h_2d->Fill(bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    //g_2d->SetPoint(g_2d->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    
                    if(bestvtx.electronTheta()<=intersecX_){
                        case_counts["t1me<m"]++;
                        case_h2d["t1me<m"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                        //case_g2d["t1me<m"]->SetPoint(case_g2d["t1me<m"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());
                    }

                    case_h1d_aco["t1mem"]->Fill(abs(bestvtx.modifiedAcoplanarity()));
                    case_h1d_vertex["t1mem"]->Fill(bestvtx.zPositionFit());
                    case_h1d_Vtxchi2["t1mem"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1mem"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1mem"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1mem"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1mem"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1mem"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1mem"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1mem");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1mem");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1mem");

                    case_h1d_aco["t1all"]->Fill(abs(bestvtx.modifiedAcoplanarity()));
                    case_h1d_vertex["t1all"]->Fill(bestvtx.zPositionFit());
                    case_h1d_Vtxchi2["t1all"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1all"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1all"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1all"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1all");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1all");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1all");

                    continue; // so skip the rest part using useTightTrackCutTgt2_ Ntrack==2
                }

                if(TGT2 && useTightTrackCutTgt2_){

                    
                    //if(abs(bestvtx.modifiedAcoplanarity()<0.3) ){
                    if( (abs(bestvtx.modifiedAcoplanarity()<0.4) && acocut) || (!acocut)  ){
                        h_2d_bstvtx->Fill(bestvtx.electronTheta(),bestvtx.muonTheta()); 
                        //g_2d_bstvtx->SetPoint(g_2d_bstvtx->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta()); 
                    }
                }
                
                for(int j=0; j<tracks.size();j++)
                {
                    if(TGT1)continue;

                    if(tracks.at(j).sector()==1) {
                        TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); in.push_back(v);muone_in.push_back(&tracks.at(j));
                        //Eugenia's cut https://indico.cern.ch/event/1476217/contributions/6217032/attachments/2962101/5210167/tesi_phd_weekly.pdf
                        //if(v.Theta()>4e-3)continue;
                        sec1++;
                    }
                    if(tracks.at(j).sector()==2 && MF) {
                        if(tracks.at(j).isMuon()){sec2muon++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); outmuon.push_back(v);muone_outmuon.push_back(&tracks.at(j));}
                        else {sec2e++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); oute.push_back(v);muone_oute.push_back(&tracks.at(j));}
                    }    
                    else if(tracks.at(j).sector()==2){
                        sec2++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out.push_back(v);muone_out.push_back(&tracks.at(j));
                    }

                }

                for(int j=0; j<tracks.size();j++)
                {
                    if(TGT2)continue;
                    
                    if(tracks.at(j).sector()==0) {
                        TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); in.push_back(v);muone_in.push_back(&tracks.at(j));
                        sec0++;
                    }
                    if(tracks.at(j).sector()==1 && MF) {
                        if(tracks.at(j).isMuon()){sec1muon++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); outmuon.push_back(v);muone_outmuon.push_back(&tracks.at(j));}
                        else {sec1e++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); oute.push_back(v);muone_oute.push_back(&tracks.at(j));}
                    }    
                    else if(tracks.at(j).sector()==1){
                        sec1++; TVector3 v(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); v=v.Unit(); out.push_back(v);muone_out.push_back(&tracks.at(j));
                    }

                }

                double angle_e, angle_mu;
                double angle0, angle1;
                double aco = 1000;

                double tmp_h1d_x[3],tmp_h1d_y[3],tmp_h1d_r[3],tmp_h1d_dx,tmp_h1d_dy,tmp_h1d_dr;

                //case1: with MF
                if (sec1 == 1 && sec2e == 1 && sec2muon == 1 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(outmuon.at(0)); 
                    aco=acoplanarity(in.at(0),oute.at(0),outmuon.at(0)); 
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
                    //if( abs(aco)>0.3)continue;//0.3 rad

                    // follow mem
                    tmp_h1d_x[0] = CalculateXtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_x[0]["t1all"]->Fill(tmp_h1d_x[0]);
                    case_h1d_x[0]["t1mem"]->Fill(tmp_h1d_x[0]);
                    tmp_h1d_x[1] = CalculateXtgt(muone_oute.at(0),z_tgt2_);
                    case_h1d_x[1]["t1all"]->Fill(tmp_h1d_x[1]);
                    case_h1d_x[1]["t1mem"]->Fill(tmp_h1d_x[1]);
                    tmp_h1d_x[2] = CalculateXtgt(muone_outmuon.at(0),z_tgt2_);
                    case_h1d_x[2]["t1all"]->Fill(tmp_h1d_x[2]);
                    case_h1d_x[2]["t1mem"]->Fill(tmp_h1d_x[2]);

                    case_h1d_bstvtx_x["t1all"]->Fill(bestvtx.xPositionFit());
                    case_h1d_bstvtx_x["t1mem"]->Fill(bestvtx.xPositionFit());

                    tmp_h1d_y[0] = CalculateYtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_y[0]["t1all"]->Fill(tmp_h1d_y[0]);
                    case_h1d_y[0]["t1mem"]->Fill(tmp_h1d_y[0]);
                    tmp_h1d_y[1] = CalculateYtgt(muone_oute.at(0),z_tgt2_);
                    case_h1d_y[1]["t1all"]->Fill(tmp_h1d_y[1]);
                    case_h1d_y[1]["t1mem"]->Fill(tmp_h1d_y[1]);
                    tmp_h1d_y[2] = CalculateYtgt(muone_outmuon.at(0),z_tgt2_);
                    case_h1d_y[2]["t1all"]->Fill(tmp_h1d_y[2]);
                    case_h1d_y[2]["t1mem"]->Fill(tmp_h1d_y[2]);

                    case_h1d_bstvtx_y["t1all"]->Fill(bestvtx.yPositionFit());
                    case_h1d_bstvtx_y["t1mem"]->Fill(bestvtx.yPositionFit());

                    tmp_h1d_r[0] = CalculateRtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_r[0]["t1all"]->Fill(tmp_h1d_r[0]);
                    case_h1d_r[0]["t1mem"]->Fill(tmp_h1d_r[0]);
                    tmp_h1d_r[1] = CalculateRtgt(muone_oute.at(0),z_tgt2_);
                    case_h1d_r[1]["t1all"]->Fill(tmp_h1d_r[1]);
                    case_h1d_r[1]["t1mem"]->Fill(tmp_h1d_r[1]);
                    tmp_h1d_r[2] = CalculateRtgt(muone_outmuon.at(0),z_tgt2_);
                    case_h1d_r[2]["t1all"]->Fill(tmp_h1d_r[2]);
                    case_h1d_r[2]["t1mem"]->Fill(tmp_h1d_r[2]);

                    case_h1d_bstvtx_r["t1all"]->Fill( sqrt( bestvtx.yPositionFit()*bestvtx.yPositionFit()+bestvtx.xPositionFit()*bestvtx.xPositionFit()) );
                    case_h1d_bstvtx_r["t1mem"]->Fill( sqrt( bestvtx.yPositionFit()*bestvtx.yPositionFit()+bestvtx.xPositionFit()*bestvtx.xPositionFit()) );

                    
                    tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[1];
                    case_h1d_dx[0]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[0]["t1mem"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[1];
                    case_h1d_dy[0]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[0]["t1mem"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[1];
                    case_h1d_dr[0]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[0]["t1mem"]->Fill(tmp_h1d_dr);

                    tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[2];
                    case_h1d_dx[1]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[1]["t1mem"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[2];
                    case_h1d_dy[1]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[1]["t1mem"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[2];
                    case_h1d_dr[1]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[1]["t1mem"]->Fill(tmp_h1d_dr);

                    tmp_h1d_dx = tmp_h1d_x[1]-tmp_h1d_x[2];
                    case_h1d_dx[2]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[2]["t1mem"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[1]-tmp_h1d_y[2];
                    case_h1d_dy[2]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[2]["t1mem"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[1]-tmp_h1d_r[2];
                    case_h1d_dr[2]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[2]["t1mem"]->Fill(tmp_h1d_dr);
                    
                    
                    case_counts["t1all"]++;
                    case_counts["t1mem"]++;
                    case_h2d["t1all"]->Fill(angle_e,angle_mu);
                    case_h2d["t1mem"]->Fill(angle_e,angle_mu);
                    case_h2d_bstvtx["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d_bstvtx["t1mem"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_e,angle_mu);
                    //case_g2d["t1mem"]->SetPoint(case_g2d["t1mem"]->GetN(),angle_e,angle_mu);
                    //case_g2d_bstvtx["t1all"]->SetPoint(case_g2d_bstvtx["t1all"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());
                    //case_g2d_bstvtx["t1mem"]->SetPoint(case_g2d_bstvtx["t1mem"]->GetN(),bestvtx.electronTheta(),bestvtx.muonTheta());

                    h_2d->Fill(angle_e,angle_mu); 
                    //g_2d->SetPoint(g_2d->GetN(),angle_e,angle_mu); 

                    if(angle_e<=intersecX_){

                        case_counts["t1me<m"]++;
                        
                        case_h2d["t1me<m"]->Fill(angle_e,angle_mu);
                        case_h2d_bstvtx["t1me<m"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                        
                        case_h1d_vertex["t1me<m"]->Fill(bestvtx.zPositionFit());

                        //case_g2d["t1me<m"]->SetPoint(case_g2d["t1me<m"]->GetN(),angle_e,angle_mu);
                        //case_g2d_bstvtx["t1me<m"]->SetPoint(case_g2d_bstvtx["t1me<m"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());
                    }
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                    //case_h2d_bstvtx["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    case_h1d_aco["t1mem"]->Fill(abs(aco));
                    case_h1d_vertex["t1mem"]->Fill(bestvtx.zPositionFit());
                    case_h1d_Vtxchi2["t1mem"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1mem"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1mem"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1mem"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1mem"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1mem"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1mem"]->Fill(LeftOverHits2.size());

                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1mem");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1mem");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1mem");

                    case_h1d_aco["t1all"]->Fill(abs(aco));
                    case_h1d_vertex["t1all"]->Fill(bestvtx.zPositionFit());
                    case_h1d_Vtxchi2["t1all"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1all"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1all"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1all"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1all");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1all");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1all");


                }    
                if (sec1 == 1 && sec2e == 2 && sec2muon == 0 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(oute.at(1));    
                    aco=acoplanarity(in.at(0),oute.at(0),oute.at(1)); 
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
                    //if( abs(aco)>0.3)continue;//0.3 rad

                    // follow mee
                    tmp_h1d_x[0] = CalculateXtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_x[0]["t1all"]->Fill(tmp_h1d_x[0]);
                    case_h1d_x[0]["t1mee"]->Fill(tmp_h1d_x[0]);
                    tmp_h1d_x[1] = CalculateXtgt(muone_oute.at(0),z_tgt2_);
                    case_h1d_x[1]["t1all"]->Fill(tmp_h1d_x[1]);
                    case_h1d_x[1]["t1mee"]->Fill(tmp_h1d_x[1]);
                    tmp_h1d_x[2] = CalculateXtgt(muone_oute.at(1),z_tgt2_);
                    case_h1d_x[2]["t1all"]->Fill(tmp_h1d_x[2]);
                    case_h1d_x[2]["t1mee"]->Fill(tmp_h1d_x[2]);

                    case_h1d_bstvtx_x["t1all"]->Fill(bestvtx.xPositionFit());
                    case_h1d_bstvtx_x["t1mee"]->Fill(bestvtx.xPositionFit());

                    tmp_h1d_y[0] = CalculateYtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_y[0]["t1all"]->Fill(tmp_h1d_y[0]);
                    case_h1d_y[0]["t1mee"]->Fill(tmp_h1d_y[0]);
                    tmp_h1d_y[1] = CalculateYtgt(muone_oute.at(0),z_tgt2_);
                    case_h1d_y[1]["t1all"]->Fill(tmp_h1d_y[1]);
                    case_h1d_y[1]["t1mee"]->Fill(tmp_h1d_y[1]);
                    tmp_h1d_y[2] = CalculateYtgt(muone_oute.at(1),z_tgt2_);
                    case_h1d_y[2]["t1all"]->Fill(tmp_h1d_y[2]);
                    case_h1d_y[2]["t1mee"]->Fill(tmp_h1d_y[2]);

                    case_h1d_bstvtx_y["t1all"]->Fill(bestvtx.yPositionFit());
                    case_h1d_bstvtx_y["t1mee"]->Fill(bestvtx.yPositionFit());

                    tmp_h1d_r[0] = CalculateRtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_r[0]["t1all"]->Fill(tmp_h1d_r[0]);
                    case_h1d_r[0]["t1mee"]->Fill(tmp_h1d_r[0]);
                    tmp_h1d_r[1] = CalculateRtgt(muone_oute.at(0),z_tgt2_);
                    case_h1d_r[1]["t1all"]->Fill(tmp_h1d_r[1]);
                    case_h1d_r[1]["t1mee"]->Fill(tmp_h1d_r[1]);
                    tmp_h1d_r[2] = CalculateRtgt(muone_oute.at(1),z_tgt2_);
                    case_h1d_r[2]["t1all"]->Fill(tmp_h1d_r[2]);
                    case_h1d_r[2]["t1mee"]->Fill(tmp_h1d_r[2]);

                    case_h1d_bstvtx_r["t1all"]->Fill( sqrt( bestvtx.yPositionFit()*bestvtx.yPositionFit()+bestvtx.xPositionFit()*bestvtx.xPositionFit()) );
                    case_h1d_bstvtx_r["t1mee"]->Fill( sqrt( bestvtx.yPositionFit()*bestvtx.yPositionFit()+bestvtx.xPositionFit()*bestvtx.xPositionFit()) );

                    
                    tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[1];
                    case_h1d_dx[0]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[0]["t1mee"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[1];
                    case_h1d_dy[0]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[0]["t1mee"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[1];
                    case_h1d_dr[0]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[0]["t1mee"]->Fill(tmp_h1d_dr);

                    tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[2];
                    case_h1d_dx[1]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[1]["t1mee"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[2];
                    case_h1d_dy[1]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[1]["t1mee"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[2];
                    case_h1d_dr[1]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[1]["t1mee"]->Fill(tmp_h1d_dr);

                    tmp_h1d_dx = tmp_h1d_x[1]-tmp_h1d_x[2];
                    case_h1d_dx[2]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[2]["t1mee"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[1]-tmp_h1d_y[2];
                    case_h1d_dy[2]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[2]["t1mee"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[1]-tmp_h1d_r[2];
                    case_h1d_dr[2]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[2]["t1mee"]->Fill(tmp_h1d_dr);
                    
                    case_counts["t1mee"]++; 
                    case_counts["t1all"]++; 

                    if(angle_e>angle_mu){
                        case_h2d["t1mee"]->Fill(angle_e,angle_mu);case_h2d["t1all"]->Fill(angle_e,angle_mu);
                        //case_g2d["t1mee"]->SetPoint(case_g2d["t1mee"]->GetN(),angle_e,angle_mu);
                        //case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t1mee"]->Fill(angle_mu,angle_e);case_h2d["t1all"]->Fill(angle_mu,angle_e);
                        //case_g2d["t1mee"]->SetPoint(case_g2d["t1mee"]->GetN(), angle_mu,angle_e);
                        //case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(), angle_mu,angle_e);
                    }
                    
                    case_h2d_bstvtx["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d_bstvtx["t1mee"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    case_h1d_aco["t1mee"]->Fill(abs(aco));
                    case_h1d_vertex["t1mee"]->Fill(bestvtx.zPositionFit());
                    case_h1d_Vtxchi2["t1mee"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1mee"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1mee"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1mee"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1mee"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1mee"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1mee"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1mee");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1mee");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1mee");


                    case_h1d_aco["t1all"]->Fill(abs(aco));
                    case_h1d_vertex["t1all"]->Fill(bestvtx.zPositionFit());
                    case_h1d_Vtxchi2["t1all"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1all"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1all"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1all"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1all");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1all");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1all");

                    

                    //case_g2d_bstvtx["t1all"]->SetPoint(case_g2d_bstvtx["t1all"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());
                    //case_g2d_bstvtx["t1mee"]->SetPoint(case_g2d_bstvtx["t1mee"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                    //case_h2d_bstvtx["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                }    
                int mmm_matching = 1;

                if (sec1 == 1 && sec2e == 0 && sec2muon == 2 && MF){ 
                    angle_e=in.at(0).Angle(outmuon.at(0)); 
                    angle_mu=in.at(0).Angle(outmuon.at(1)); 
                    aco=acoplanarity(in.at(0),outmuon.at(0),outmuon.at(1)); 
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
                    //if( abs(aco)>0.3)continue;//0.3 rad
                    
                    // follow mmm
                    tmp_h1d_x[0] = CalculateXtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_x[0]["t1all"]->Fill(tmp_h1d_x[0]);
                    case_h1d_x[0]["t1mmm"]->Fill(tmp_h1d_x[0]);
                    
                    tmp_h1d_x[1] = CalculateXtgt(muone_outmuon.at(0),z_tgt2_);
                    tmp_h1d_x[2] = CalculateXtgt(muone_outmuon.at(1),z_tgt2_);

                    case_h1d_x[1]["t1all"]->Fill(tmp_h1d_x[1]);
                    case_h1d_x[1]["t1mmm"]->Fill(tmp_h1d_x[1]);
                    case_h1d_x[2]["t1all"]->Fill(tmp_h1d_x[2]);
                    case_h1d_x[2]["t1mmm"]->Fill(tmp_h1d_x[2]);

                    case_h1d_bstvtx_x["t1all"]->Fill(bestvtx.xPositionFit());
                    case_h1d_bstvtx_x["t1mmm"]->Fill(bestvtx.xPositionFit());



                    tmp_h1d_y[0] = CalculateYtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_y[0]["t1all"]->Fill(tmp_h1d_y[0]);
                    case_h1d_y[0]["t1mmm"]->Fill(tmp_h1d_y[0]);
                    
                    tmp_h1d_y[1] = CalculateYtgt(muone_outmuon.at(0),z_tgt2_);
                    tmp_h1d_y[2] = CalculateYtgt(muone_outmuon.at(1),z_tgt2_);

                    
                    case_h1d_y[1]["t1all"]->Fill(tmp_h1d_y[1]);
                    case_h1d_y[1]["t1mmm"]->Fill(tmp_h1d_y[1]);
                    
                    case_h1d_y[2]["t1all"]->Fill(tmp_h1d_y[2]);
                    case_h1d_y[2]["t1mmm"]->Fill(tmp_h1d_y[2]);

                    case_h1d_bstvtx_y["t1all"]->Fill(bestvtx.yPositionFit());
                    case_h1d_bstvtx_y["t1mmm"]->Fill(bestvtx.yPositionFit());

                    tmp_h1d_r[0] = CalculateRtgt(muone_in.at(0),z_tgt2_);
                    case_h1d_r[0]["t1all"]->Fill(tmp_h1d_r[0]);
                    case_h1d_r[0]["t1mmm"]->Fill(tmp_h1d_r[0]);
                    
                    tmp_h1d_r[1] = CalculateRtgt(muone_outmuon.at(0),z_tgt2_);
                    tmp_h1d_r[2] = CalculateRtgt(muone_outmuon.at(1),z_tgt2_);

                    case_h1d_r[1]["t1all"]->Fill(tmp_h1d_r[1]);
                    case_h1d_r[1]["t1mmm"]->Fill(tmp_h1d_r[1]);
                    case_h1d_r[2]["t1all"]->Fill(tmp_h1d_r[2]);
                    case_h1d_r[2]["t1mmm"]->Fill(tmp_h1d_r[2]);

                    case_h1d_bstvtx_r["t1all"]->Fill( sqrt( bestvtx.yPositionFit()*bestvtx.yPositionFit()+bestvtx.xPositionFit()*bestvtx.xPositionFit()) );
                    case_h1d_bstvtx_r["t1mmm"]->Fill( sqrt( bestvtx.yPositionFit()*bestvtx.yPositionFit()+bestvtx.xPositionFit()*bestvtx.xPositionFit()) );

                    
                    tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[1];
                    if(mmm_matching &&  abs(tmp_h1d_dx)> abs(tmp_h1d_x[0]-tmp_h1d_x[2]) )tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[2];
                    case_h1d_dx[0]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[0]["t1mmm"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[1];
                    if(mmm_matching &&  abs(tmp_h1d_dy)> abs(tmp_h1d_y[0]-tmp_h1d_y[2]) )tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[2];
                    case_h1d_dy[0]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[0]["t1mmm"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[1];
                    if(mmm_matching &&  abs(tmp_h1d_dr)> abs(tmp_h1d_r[0]-tmp_h1d_r[2]) )tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[2];
                    case_h1d_dr[0]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[0]["t1mmm"]->Fill(tmp_h1d_dr);

                    tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[2];
                    if(mmm_matching &&  abs(tmp_h1d_dx)< abs(tmp_h1d_x[0]-tmp_h1d_x[1]) )tmp_h1d_dx = tmp_h1d_x[0]-tmp_h1d_x[1];
                    case_h1d_dx[1]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[1]["t1mmm"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[2];
                    if(mmm_matching &&  abs(tmp_h1d_dy)< abs(tmp_h1d_y[0]-tmp_h1d_y[1]) )tmp_h1d_dy = tmp_h1d_y[0]-tmp_h1d_y[1];
                    case_h1d_dy[1]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[1]["t1mmm"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[2];
                    if(mmm_matching &&  abs(tmp_h1d_dr)< abs(tmp_h1d_r[0]-tmp_h1d_r[1]) )tmp_h1d_dr = tmp_h1d_r[0]-tmp_h1d_r[1];
                    case_h1d_dr[1]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[1]["t1mmm"]->Fill(tmp_h1d_dr);

                    tmp_h1d_dx = tmp_h1d_x[1]-tmp_h1d_x[2];
                    case_h1d_dx[2]["t1all"]->Fill(tmp_h1d_dx);
                    case_h1d_dx[2]["t1mmm"]->Fill(tmp_h1d_dx);
                    tmp_h1d_dy = tmp_h1d_y[1]-tmp_h1d_y[2];
                    case_h1d_dy[2]["t1all"]->Fill(tmp_h1d_dy);
                    case_h1d_dy[2]["t1mmm"]->Fill(tmp_h1d_dy);
                    tmp_h1d_dr = tmp_h1d_r[1]-tmp_h1d_r[2];
                    case_h1d_dr[2]["t1all"]->Fill(tmp_h1d_dr);
                    case_h1d_dr[2]["t1mmm"]->Fill(tmp_h1d_dr);

                    case_counts["t1mmm"]++; 
                    case_counts["t1all"]++; 

                    
                    //if( bestvtx.muonTheta()>0.0003 && (bestvtx.zPositionFit()<600 || bestvtx.zPositionFit()>1200))std::cout<<"checking strange mmm: "<<bestvtx.zPositionFit()<<std::endl;
                    
                    if(angle_e>angle_mu){
                        case_h2d["t1mmm"]->Fill(angle_e,angle_mu);case_h2d["t1all"]->Fill(angle_e,angle_mu);
                        //case_g2d["t1mmm"]->SetPoint(case_g2d["t1mmm"]->GetN(),angle_e,angle_mu);
                        //case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t1mmm"]->Fill(angle_mu,angle_e);case_h2d["t1all"]->Fill(angle_mu,angle_e);
                        //case_g2d["t1mmm"]->SetPoint(case_g2d["t1mmm"]->GetN(),angle_mu,angle_e);
                        //case_g2d["t1all"]->SetPoint(case_g2d["t1all"]->GetN(),angle_mu,angle_e);
                    }
                    
                    case_h2d_bstvtx["t1all"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    case_h2d_bstvtx["t1mmm"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                    
                    //case_g2d_bstvtx["t1all"]->SetPoint(case_g2d_bstvtx["t1all"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());
                    //case_g2d_bstvtx["t1mmm"]->SetPoint(case_g2d_bstvtx["t1mmm"]->GetN(), bestvtx.electronTheta(),bestvtx.muonTheta());

                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                    //case_h2d_bstvtx["golden"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());

                    case_h1d_aco["t1mmm"]->Fill(abs(aco));
                    if(bestvtx.zPositionFit()>0 && bestvtx.zPositionFit()<1200)case_h1d_vertex["t1mmm"]->Fill(bestvtx.zPositionFit());
                    else if(bestvtx.zPositionFit()<0)case_h1d_vertex["t1mmm"]->Fill(1);
                    else if(bestvtx.zPositionFit()>1200)case_h1d_vertex["t1mmm"]->Fill(1119);
                    case_h1d_Vtxchi2["t1mmm"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1mmm"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1mmm"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1mmm"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1mmm"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1mmm"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1mmm"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1mmm");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1mmm");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1mmm");

                    case_h1d_aco["t1all"]->Fill(abs(aco));
                    if(bestvtx.zPositionFit()>0 && bestvtx.zPositionFit()<1200)case_h1d_vertex["t1all"]->Fill(bestvtx.zPositionFit());
                    else if(bestvtx.zPositionFit()<0)case_h1d_vertex["t1all"]->Fill(1);
                    else if(bestvtx.zPositionFit()>1200)case_h1d_vertex["t1all"]->Fill(1119);
                    case_h1d_Vtxchi2["t1all"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                    case_h1d_LeftOverhits0["t1all"]->Fill(LeftOverHits0.size());
                    case_h1d_LeftOverhits1["t1all"]->Fill(LeftOverHits1.size());
                    case_h1d_LeftOverhits2["t1all"]->Fill(LeftOverHits2.size());
                    //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits0.size());
                    //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits1.size());
                    //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1all"]->Fill(LeftOverHits2.size());
                    processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1all");
                    processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1all");
                    processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1all");

                    if( bestvtx.muonTheta()>0.0003 ){
                        case_h2d_bstvtx["t1mmmband"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                        if(bestvtx.zPositionFit()>0 && bestvtx.zPositionFit()<1200)case_h1d_vertex["t1mmmband"]->Fill(bestvtx.zPositionFit());
                        else if(bestvtx.zPositionFit()<0)case_h1d_vertex["t1mmmband"]->Fill(1);
                        else if(bestvtx.zPositionFit()>1200)case_h1d_vertex["t1mmmband"]->Fill(1119);
                        case_h1d_aco["t1mmmband"]->Fill(abs(aco));
                        case_h1d_Vtxchi2["t1mmmband"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                        case_h1d_LeftOverhits0["t1mmmband"]->Fill(LeftOverHits0.size());
                        case_h1d_LeftOverhits1["t1mmmband"]->Fill(LeftOverHits1.size());
                        case_h1d_LeftOverhits2["t1mmmband"]->Fill(LeftOverHits2.size());
                        //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1mmmband"]->Fill(LeftOverHits0.size());
                        //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1mmmband"]->Fill(LeftOverHits1.size());
                        //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1mmmband"]->Fill(LeftOverHits2.size());
                        processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1mmmband");
                        processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1mmmband");
                        processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1mmmband");

                        if(angle_e>angle_mu){
                            case_h2d["t1mmmband"]->Fill(angle_e,angle_mu);
                        }
                        else {
                            case_h2d["t1mmmband"]->Fill(angle_mu,angle_e);
                        }
                        

                    }
                    else{

                        case_h2d_bstvtx["t1mmmoutofband"]->Fill(bestvtx.electronTheta(),bestvtx.muonTheta());
                        if(bestvtx.zPositionFit()>0 && bestvtx.zPositionFit()<1200)case_h1d_vertex["t1mmmoutofband"]->Fill(bestvtx.zPositionFit());
                        else if(bestvtx.zPositionFit()<0)case_h1d_vertex["t1mmmoutofband"]->Fill(1);
                        else if(bestvtx.zPositionFit()>1200)case_h1d_vertex["t1mmmoutofband"]->Fill(1119);
                        case_h1d_aco["t1mmmoutofband"]->Fill(abs(aco));
                        case_h1d_Vtxchi2["t1mmmoutofband"]->Fill(bestvtx.chi2perDegreeOfFreedom());
                        case_h1d_LeftOverhits0["t1mmmoutofband"]->Fill(LeftOverHits0.size());
                        case_h1d_LeftOverhits1["t1mmmoutofband"]->Fill(LeftOverHits1.size());
                        case_h1d_LeftOverhits2["t1mmmoutofband"]->Fill(LeftOverHits2.size());
                        //for (auto const& hit : LeftOverHits0)case_h1d_LeftOverhits0perModule[hit->moduleID()]["t1mmmoutofband"]->Fill(LeftOverHits0.size());
                        //for (auto const& hit : LeftOverHits1)case_h1d_LeftOverhits1perModule[hit->moduleID()]["t1mmmoutofband"]->Fill(LeftOverHits1.size());
                        //for (auto const& hit : LeftOverHits2)case_h1d_LeftOverhits2perModule[hit->moduleID()]["t1mmmoutofband"]->Fill(LeftOverHits2.size());
                        processLeftoverHits(LeftOverHits0, case_h1d_LeftOverhits0perModule, "t1mmmoutofband");
                        processLeftoverHits(LeftOverHits1, case_h1d_LeftOverhits1perModule, "t1mmmoutofband");
                        processLeftoverHits(LeftOverHits2, case_h1d_LeftOverhits2perModule, "t1mmmoutofband");

                        if(angle_e>angle_mu){
                            case_h2d["t1mmmoutofband"]->Fill(angle_e,angle_mu);
                        }
                        else {
                            case_h2d["t1mmmoutofband"]->Fill(angle_mu,angle_e);
                        }

                    }

                }    
                if (sec0 == 1 && sec1e == 1 && sec1muon == 1 && MF){ 

                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(outmuon.at(0)); 
                    aco=acoplanarity(in.at(0),oute.at(0),outmuon.at(0)); 
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
                    
                    case_counts["t0mem"]++; 
                    case_counts["t0all"]++; 
                    case_h2d["t0mem"]->Fill(angle_e,angle_mu); 
                    case_h2d["t0all"]->Fill(angle_e,angle_mu); 
                    //case_g2d["t0mem"]->SetPoint(case_g2d["t0mem"]->GetN(), angle_e,angle_mu); 
                    //case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_e,angle_mu); 
                    
                    case_h1d_vertex["t0mem"]->Fill(bestvtx.zPositionFit());
                    case_h1d_vertex["t0all"]->Fill(bestvtx.zPositionFit());

                    h_2d->Fill(angle_e,angle_mu); 
                    //g_2d->SetPoint(g_2d->GetN(), angle_e,angle_mu); 
                    
                    if(angle_e<=intersecX_){
                        case_counts["t0me<m"]++;
                        case_h2d["t0me<m"]->Fill(angle_e,angle_mu); 
                        case_h1d_vertex["t0me<m"]->Fill(bestvtx.zPositionFit());
                    
                        //case_g2d["t0me<m"]->SetPoint(case_g2d["t0me<m"]->GetN(),angle_e,angle_mu); 
                    }
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                }    
                if (sec0 == 1 && sec1e == 2 && sec1muon == 0 && MF){ 
                    angle_e=in.at(0).Angle(oute.at(0));    
                    angle_mu=in.at(0).Angle(oute.at(1));    
                    aco=acoplanarity(in.at(0),oute.at(0),oute.at(1)); 
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
                    
                    case_counts["t0mee"]++; 
                    case_counts["t0all"]++; 

                    case_h1d_vertex["t0mee"]->Fill(bestvtx.zPositionFit());
                    case_h1d_vertex["t0all"]->Fill(bestvtx.zPositionFit());

                    if(angle_e>angle_mu){
                        case_h2d["t0mee"]->Fill(angle_e,angle_mu);case_h2d["t0all"]->Fill(angle_e,angle_mu);
                        //case_g2d["t0mee"]->SetPoint(case_g2d["t0mee"]->GetN(), angle_e,angle_mu);
                        //case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t0mee"]->Fill(angle_mu,angle_e);case_h2d["t0all"]->Fill(angle_mu,angle_e);
                        //case_g2d["t0mee"]->SetPoint(case_g2d["t0mee"]->GetN(), angle_mu,angle_e);
                        //case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_mu,angle_e);
                    }
                    //case_h2d["golden"]->Fill(angle_e,angle_mu);
                }    
                if (sec0 == 1 && sec1e == 0 && sec1muon == 2 && MF){ 
                    angle_e=in.at(0).Angle(outmuon.at(0)); 
                    angle_mu=in.at(0).Angle(outmuon.at(1)); 
                    aco=acoplanarity(in.at(0),outmuon.at(0),outmuon.at(1)); 
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
                    
                    case_counts["t0mmm"]++; 
                    case_counts["t0all"]++; 

                    case_h1d_vertex["t0mmm"]->Fill(bestvtx.zPositionFit());
                    case_h1d_vertex["t0all"]->Fill(bestvtx.zPositionFit());
                    
                    if(angle_e>angle_mu){
                        case_h2d["t0mmm"]->Fill(angle_e,angle_mu);case_h2d["t0all"]->Fill(angle_e,angle_mu);
                        //case_g2d["t0mmm"]->SetPoint(case_g2d["t0mmm"]->GetN(), angle_e,angle_mu);
                        //case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_e,angle_mu);
                    }
                    else {
                        case_h2d["t0mmm"]->Fill(angle_mu,angle_e);case_h2d["t0all"]->Fill(angle_mu,angle_e);
                        //case_g2d["t0mmm"]->SetPoint(case_g2d["t0mmm"]->GetN(), angle_mu,angle_e);
                        //case_g2d["t0all"]->SetPoint(case_g2d["t0all"]->GetN(), angle_mu,angle_e);
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
                    
                    if( abs(aco)>0.4 && acocut)continue;//0.4 rad
            
                    //if(tracks.size()!=3 || (angle0>angle1 && angle0>0.032) || (angle0<angle1 && angle1>0.032) || (angle0>angle1 && angle1<0.0002) || (angle0<angle1 && angle0<0.0002) )continue;
                    //flag_good_event = 1;

                    if(angle0>angle1) {
                        h_2d->Fill(angle0,angle1);
                        //g_2d->SetPoint(g_2d->GetN(),angle0,angle1);
                    }
                    else {
                        h_2d->Fill(angle1,angle0);
                        //g_2d->SetPoint(g_2d->GetN(),angle1,angle0);
                    }
                    
                }

            }//isGolden

        }// simple N tracks cut
    
    }//event loop

    std::cout<<"goldenevents: "<<goldenevents_<<"/"<<N<<std::endl;
}
