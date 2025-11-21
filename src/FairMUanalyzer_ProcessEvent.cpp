#include "FairMUanalyzer.h"

void FairMUanalyzer::ProcessEvent(Long64_t i, int debug) {
    
    int TGT1 = tgt_==0?1:0;
    int TGT2 = tgt_==1?1:0;
    //std::cout<<"tgt_"<<tgt_<<std::endl;
    if(!TGT1 && !TGT2){std::cout<<"wrong set up on TGT! return."<<std::endl;return;}

    isGolden_ = false;  

    cbmsim_->GetEntry(i);
    cbmsim_->SetBranchAddress("ReconstructionOutput", &reco_);

    const auto& tracks = reco_->reconstructedTracks();
    const auto& hits = reco_->reconstructedHits();
    const auto& bestvtx = reco_->bestVertex();

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

    bool isGolden = true;

    if(debug){
        std::cout<<"tracks.size(): "<<tracks.size()<<std::endl;
        int ntrk_sec0=0; 
        int ntrk_sec1=0;
        int ntrk_sec2=0;
        for (auto const& track : tracks) {
            std::set<int> modules;
            for (auto const& h : track.hits()) {
                modules.insert(h.moduleID());
            }
            if(debug)std::cout<<"modules.size(): "<<modules.size()<<std::endl;

            if(track.sector()==0)ntrk_sec0++;
            if(track.sector()==1)ntrk_sec1++;
            if(track.sector()==2)ntrk_sec2++;
        }
        
        std::cout<<"ntrk_sec012: "<<ntrk_sec0<<ntrk_sec1<<ntrk_sec2<<std::endl;
    }

    if ( (TGT2 && tracks.size() == 4) || (TGT1 && tracks.size() >= 3) ) {

        std::set<int> sectors;
        std::set<int> sectors01;

        int ntrk_sec0=0; 
        int ntrk_sec1=0;
        int ntrk_sec2=0;

        for (auto const& track : tracks) {
            std::set<int> modules;
            
            // need a further dealing - 21Nov25 - see chat GPT msg, for new 0.17.6 version
            
            //for (auto const& h : track.hits()) {
            //    modules.insert(h.moduleID());
            //}
            
            if(debug)std::cout<<"modules.size(): "<<modules.size()<<std::endl;
            if (modules.size() != 6) {
                //golden muon step #2: 1 hit/module
                isGolden = false;
                break;
            }
            


            //golden muon step #3: reduced chi2
            // changed 21Nov2025 - not necessary 
            //if(track.chi2perDegreeOfFreedom()>=2)isGolden = false;

            sectors.insert(track.sector());
            if(track.sector()<2)sectors01.insert(track.sector());

            if(track.sector()==0)ntrk_sec0++;
            if(track.sector()==1)ntrk_sec1++;
            if(track.sector()==2)ntrk_sec2++;

            

        }

        if(debug){
            std::cout<<"ntrk_sec012: "<<ntrk_sec0<<ntrk_sec1<<ntrk_sec2<<std::endl;
        }

        //golden muon step #4: 1/2 tracks/station, and nhit_giovanni cut
        //note we don't require muonid at this moment, just (1track + 2track + (1+) tracks) or (1track + 1track + 2 tracks) signature is ok
        //we have to do this way because we  need those 'wrong' events (even if we know by MF) to for the next step to see the count/distributions



        if(TGT1
            && 
            ( sectors01.size() != 2 || !(ntrk_sec0==1 && ntrk_sec1==2) || (nhits_sec1 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
            )isGolden = false;
            

        if(TGT2 
            && 
            ( sectors.size() != 3 || !(ntrk_sec0==1 && ntrk_sec1==1 && ntrk_sec2==2) || (nhits_sec2 > maxNhitInStat_ ) ) // suggested by Giovanni A, can be tested by turning HitCutsOn ON/OFF
            )isGolden = false;
    }
    else isGolden = false;

    isGolden_ = isGolden;

    int n_muons = 0;
    std::vector<const MUonERecoOutputTrackAnalysis*> muon_tracks;
    for (auto const& track : tracks) {
        if (track.isMuon() ) {
            n_muons++;
            muon_tracks.push_back(&track);
        }
    }

    //n_muons_ = n_muons;

    /*
    std::set<int> modules;

    for (const auto& track : tracks) {

        if (!track.isMuon()) continue;
        std::set<int> thisTrackModules;

        for (int j = 0; j < track.nHits(); ++j) {
            const auto& hit = track.GetHit(j);
            thisTrackModules.insert(hit.moduleID());
        }

        if (thisTrackModules.size() == 6) {
            modules.insert(thisTrackModules.begin(), thisTrackModules.end());
            n_muons++;
        }
    }

    if (n_muons == 3 && modules.size() == 6) {
        isGolden_ = true;
    }
    */
}