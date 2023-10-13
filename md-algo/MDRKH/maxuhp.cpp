#include "struct.hpp"

struct compare_for_nth{
    bool operator() (SItem const& a, SItem const& b) const{
        return a.sco > b.sco;
    }
};

map<KSet,KI> compute_all_ksets(DataSet const &SD, DataSet const &kSD, vector<WVector> const& WVS, int const& k){

    map<KSet,KI> KSS;
    for(int i = 0; i < WVS.size(); i++){
        vector<SItem> SIL(SD.size() + kSD.size());
        for(int j = 0; j < SD.size(); j++){
            SIL[j] = {j,get_score(SD[j], WVS[i])};
        }
        for(int j = 0; j < kSD.size(); j++){
            SIL[j+SD.size()] = {-1,get_score(kSD[j], WVS[i])};
        }
        nth_element(SIL.begin(), SIL.begin()+k-1, SIL.end(), compare_for_nth());

        KSet kset;
        for(int j = 0; j < k; j++) if(SIL[j].sid != -1) kset.insert(SIL[j].sid);

        if(KSS.count(kset) == 0) KSS.insert({kset, {WVS[i],1}});
        else KSS[kset].ct += 1;
    }

    return KSS;
}

void initial_data_structures(map<KSet,KI> const& KSS,
    vector<KSet> & KSV, vector<WVector> & CWV, vector<int> & KCT, vector<set<KID>> & TC, map<SID, int> &SCT){

    KID kid = 0;
    for(auto const& mpair: KSS){
        KSV[kid] = mpair.first;
        CWV[kid] = mpair.second.wv;
        KCT[kid] = mpair.second.ct;
        for(SID const& sid: KSV[kid]){
            TC[sid].insert(kid);//inversed linker
            if(SCT.count(sid) == 0) SCT.insert({sid, KCT[kid]});
            else SCT[sid] += KCT[kid];
        }
        kid += 1;
    }
}

Sol perform_greedy(int const& r, vector<KSet> const& KSV, vector<int> const& KCT, vector<set<KID>> const& TC, 
    map<SID, int> &SCT, vector<bool> & UC){
    Sol sol;
    while(sol.size() < r){
        int mcn = 0; //the maximum cover number
        SID csid = -1; //corresponding sid
        for(auto const& mpair: SCT){
            if(mpair.second > mcn){
                csid = mpair.first;
                mcn = mpair.second;
            }
        }

        if(csid == -1) break; //sol.size() < r

        for(KID const& kid: TC[csid]){
            if(UC[kid]){
                for(SID const& sid: KSV[kid]){
                    if(SCT[sid] == KCT[kid]) SCT.erase(sid);
                    else SCT[sid] -= KCT[kid];
                }
                UC[kid] = false;
            }
        }
        sol.insert(csid);
    }
    return sol;
}

Sol perform_clustering(DataSet const &SD, int const& r, Sol const& sol, int& round,
    vector<WVector> const& CWV, vector<int> const& KCT, vector<set<KID>> const& TC){
    int kn = CWV.size();
    int sn = SD.size();

    vector<SID> solv(sol.begin(),sol.end());
    if(solv.size() == r){
        bool isupdated;
        do{
            isupdated = false;
            vector<int> PWVS(kn,-1); //partition all weight vectors to r clusers
            for(KID kid = 0; kid < kn; kid++){
                Score maxs = 0;
                int ci = -1;
                for(int i = 0; i < r; i++){
                    Score sco = get_score(SD[solv[i]], CWV[kid]);
                    if(sco > maxs){
                        maxs = sco;
                        ci = i;
                    }
                }
                PWVS[kid] = ci;
            }

            vector<vector<int>> PSCT(sn, vector<int>(r, 0)); //from sid and clusterid (in {0,..,r-1}) to partial hitting pro
            for(SID sid = 0; sid < sn; sid++){
                for(KID const& kid: TC[sid]){
                    PSCT[sid][PWVS[kid]] += KCT[kid];
                    //kid belongs to the PWVS[kid]]-th cluster
                }
            }

            vector<int> maxcts(r);
            for(int i = 0; i < r; i++) maxcts[i] = PSCT[solv[i]][i];
            for(SID sid = 0; sid < sn; sid++){
                for(int i = 0; i < r; i++){
                    if(PSCT[sid][i] > maxcts[i]){
                        maxcts[i] = PSCT[sid][i];
                        solv[i] = sid;
                        isupdated = true;
                    }
                }
            }
            round++;
        }while(isupdated);
    }
    return set<SID>(solv.begin(), solv.end());
}

