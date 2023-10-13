#include "struct.hpp"

Sol set_cover(vector<set<WID>> const& STC, int const& sn, int const& sd){
    Sol sol;
    vector<set<SID>> WVC(sd); //from WID to the set of skyline tuples covered WID weight vector
    vector<bool> UC(sd, true); //mark uncovered wvectors
    map<SID, int> CT; //the counter of the number of wvectors covered by each tuple

    for(SID i = 0; i < sn; i++){
        for(WID const& wid: STC[i]) WVC[wid].insert(i);
        if(STC[i].size() > 0) CT.insert({i, STC[i].size()});
    }

    while(CT.size() > 0){
        int mcn = 0; //maximum number of covered wvectors
        SID csid = -1; //corresponding SID 
        for(auto const& mpair: CT){
            if(mpair.second > mcn){
                csid = mpair.first;
                mcn = mpair.second;
            }
        }
        for(WID const& wid: STC[csid]){
            if(UC[wid]){
                UC[wid] = false;
                for(SID const& sid: WVC[wid]){
                    if(CT[sid] == 1) CT.erase(sid);
                    else CT[sid] -= 1; 
                }
            }
        }
        sol.insert(csid);
    }
    return sol;
}

Sol mrst_oracle(vector<vector<RRatio>> const& MAT, RRatio const& rrb, int const& sn, int const& sd){
    set<set<WID>> SEEN;
    vector<set<WID>> DSTC(sn); //the weight vectors of duplicated skyline tuple covering
    for(SID i = 0; i < sn; i++){
        set<WID> wvset;
        for(WID j = 0; j < sd; j++){
            if(MAT[i][j] <= rrb) wvset.insert(j);
        }
        if(SEEN.count(wvset) == 0){
            DSTC[i] = wvset;
            SEEN.insert(wvset);
        }
    }
    return set_cover(DSTC, sn, sd);
}