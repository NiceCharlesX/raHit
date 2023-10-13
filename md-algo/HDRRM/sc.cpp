#include "struct.hpp"

Sol set_cover(KSetS const& KSS, int const& n){
    Sol sol;
    vector<KSet> KSV; //the vector of ksets
    vector<bool> UC(KSS.size(), true); //mark uncovered ksets
    vector<set<int>> TC(n); //kset ids covered by each tuple, inverse list of KSV
    map<TID, int> CT; //the counter of the number of ksets covered by each tuple

    KID kid = 0;
    for(KSet const& kset: KSS){
        KSV.push_back(kset);
        for(TID const& tid : kset){
            TC[tid].insert(kid);
            if(CT.count(tid) == 0) CT.insert({tid, 1});
            else CT[tid] += 1;
        }
        kid += 1;
    }

    while(CT.size() > 0){
        int mcn = 0; //maximum number of covered k-sets
        TID ctid = -1; //corresponding tid
        for(auto const& mpair: CT){
            if(mpair.second > mcn){
                ctid = mpair.first;
                mcn = mpair.second;
            }
        }
        for(KID const& kid: TC[ctid]){
            if(UC[kid]){ //kset is uncovered before
                UC[kid] = false;
                for(TID const& tid: KSV[kid]){ //update counter CT
                    if(CT[tid] == 1) CT.erase(tid);
                    else CT[tid] -= 1; 
                }
            }
        }
        sol.insert(ctid);
    }
    return sol;
}