#include "struct.hpp"

struct compare_for_nth{
    bool operator() (SItem const& a, SItem const& b) const{
        return a.sco > b.sco;
    }
};

//by sum in descending order
struct compare_tuples_by_sum{
    bool operator() (Tuple const& a, Tuple const& b) const{
        Coord sa = 0, sb = 0;
        for(int i = 0; i < a.size(); i++){
            sa += a[i];
            sb += b[i];
        }
        return sa > sb;
    }
};

KSetS init_sitem_lists(vector<vector<SItem>> & SILS, vector<WVector> const& WVS, DataSet const &D, 
    Sol const& basis){
    KSetS KSS;
    for(int i = 0; i < WVS.size(); i++){
        Score maxs = 0; //max score of basis
        for(TID const& tid: basis) maxs = max(maxs, get_score(D[tid],WVS[i]));
        for(int j = 0; j < D.size(); j++){
            Score sco = get_score(D[j], WVS[i]);
            if(sco > maxs) SILS[i].push_back({j,sco}); //use basis to filter
        }
        if(SILS[i].size() > 0){
            nth_element(SILS[i].begin(), SILS[i].begin(), SILS[i].end(), compare_for_nth());
            KSet kset;
            kset.insert(SILS[i][0].tid);
            if(KSS.count(kset) == 0) KSS.insert(kset);
        }
    }
    return KSS;
}

/*
 * The first k/2 tuples are already partial sorted
 */
KSetS incre_compute_ksets(vector<vector<SItem>> & SILS, Rank const& k){
    KSetS KSS;
    for(int i = 0; i < SILS.size(); i++){
        if(SILS[i].size() >= k){ //otherwise, the weight vector is covered by basis.
            KSet kset;
            nth_element(SILS[i].begin()+k/2, SILS[i].begin()+k-1, SILS[i].end(), compare_for_nth());
            for(int j = 0; j < k; j++) kset.insert(SILS[i][j].tid);
            if(KSS.count(kset) == 0) KSS.insert(kset);
        }
    }
    return KSS;
}

/* 
 * The first #(low-1) tuples are already partial sorted
 * The first #high tuples are already partial sorted
 */ 
KSetS decre_compute_ksets(vector<vector<SItem>> & SILS, Rank const& k, Rank const& low, Rank const& high){
    KSetS KSS;
    for(int i = 0; i < SILS.size(); i++){
        if(SILS[i].size() >= k){
            KSet kset;
            if(SILS[i].size() > high){
                nth_element(SILS[i].begin()+low-1, SILS[i].begin()+k-1, SILS[i].begin()+high, compare_for_nth());
            }else{
                nth_element(SILS[i].begin()+low-1, SILS[i].begin()+k-1, SILS[i].end(), compare_for_nth());
            }
            for(int j = 0; j < k; j++) kset.insert(SILS[i][j].tid);
            if(KSS.count(kset) == 0) KSS.insert(kset);
        }
    }
    return KSS;
}