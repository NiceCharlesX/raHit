#include "struct.hpp"

/*
 * determine whether a is dominated by b
 */
bool is_dominated(Tuple const& a, Tuple const & b){
    for(int i = 0; i < a.size(); i++){
        if(a[i] > b[i]) return false;
    }
    return true;
}

/* 
 * Incremental update dominated times
 * If dtimes > k, return dtimes = k + 1
 * CLD: current level dataset
 */
bool incre_update_dtimes(Tuple const& t, DataSet const& CLD, int const& k, int & dtimes){
    bool notdom = true;
    for(Tuple const& kst: CLD){
        if(is_dominated(t, kst)){
            dtimes++;
            notdom = false;
            if(dtimes >= k) break;
        }
    }
    return notdom;
}

/* 
 * Divide the kskyband (tuple with dtimes <= k-1) into k layers
 * tuples of the same layer have no dominance relationship
 * tuples of the i-th layer constitute the skyline of the dataset after removing the first i-1 layer of tuples
 */
vector<DataSet> compute_skyband(DataSet const& D, int const& k){
    vector<DataSet> LD(k); //k-level dataset, LD[0] is the Skyline
    for(Tuple const& t: D){
        int dtimes = 0; //the times of t being dominated
        for(int i = 0; i < k; i++){
            if(incre_update_dtimes(t, LD[i], k, dtimes)){
                LD[i].push_back(t);
                break;
            }
            if(dtimes >= k) break;
        }
    }
    return LD;
}