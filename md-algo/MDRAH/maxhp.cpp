#include "struct.hpp"

map<ASet,int> compute_all_asets(DataSet const &SD, vector<WVector> const& WVS, double const& alpha){
    map<ASet,int> ASC; //count of ASet
    for(WVector const& wv: WVS){
        Score MaxS = 0;
        vector<Score> SL(SD.size()); //Score list
        for(SID sid = 0; sid < SD.size(); sid++){
            SL[sid] = get_score(SD[sid], wv);
            MaxS = max(MaxS, SL[sid]);
        }

        ASet aset;
        for(SID sid = 0; sid < SD.size(); sid++){
            if(SL[sid] >= MaxS * alpha) aset.insert(sid);
        }

        //deduplication by counting
        if(ASC.count(aset) == 0) ASC.insert({aset, 1});
        else ASC[aset] += 1;
    }
    return ASC;
}

map<ASet,int> compute_all_asets(DataSet const &SD, map<WVector, int> const& WVC, double const& alpha){
    map<ASet,int> ASC; //count of ASet
    for(auto const& mpair: WVC){
        Score MaxS = 0;
        vector<Score> SL(SD.size()); //Score list
        for(SID sid = 0; sid < SD.size(); sid++){
            SL[sid] = get_score(SD[sid], mpair.first);
            MaxS = max(MaxS, SL[sid]);
        }

        ASet aset;
        for(SID sid = 0; sid < SD.size(); sid++){
            if(SL[sid] >= MaxS * alpha) aset.insert(sid);
        }

        //deduplication by counting
        if(ASC.count(aset) == 0) ASC.insert({aset, mpair.second});
        else ASC[aset] += mpair.second;
    }
    return ASC;
}

void initial_data_structures(map<ASet,int> const& ASC,
    vector<ASet> & ASV, vector<int> & ACT, vector<set<AID>> & TC, map<SID, int> &SCT){
    AID aid = 0;
    for(auto const& mpair: ASC){
        ASV[aid] = mpair.first;
        ACT[aid] = mpair.second;
        for(SID const& sid: ASV[aid]){
            TC[sid].insert(aid);//inversed linker
            if(SCT.count(sid) == 0) SCT.insert({sid, ACT[aid]});
            else SCT[sid] += ACT[aid];
        }
        aid += 1;
    }
}

Sol perform_greedy(int const& r, vector<ASet> const& ASV, vector<int> const& ACT, vector<set<AID>> const& TC, 
    map<SID, int> &SCT, vector<bool> & UC){
    Sol sol;
    while(sol.size() < r){
        int MaxCN = 0; //the maximum cover number
        SID MaxSID = -1; //corresponding sid
        for(auto const& mpair: SCT){
            if(mpair.second > MaxCN){
                MaxSID = mpair.first;
                MaxCN = mpair.second;
            }
        }

        if(MaxSID == -1) break; //sol.size() < r

        for(AID const& aid: TC[MaxSID]){
            if(UC[aid]){
                for(SID const& sid: ASV[aid]){
                    if(SCT[sid] == ACT[aid]) SCT.erase(sid);
                    else SCT[sid] -= ACT[aid];
                }
                UC[aid] = false;
            }
        }
        sol.insert(MaxSID);
    }
    return sol;
}
