#include "struct.hpp"

HP test_solahp(DataSet const& D, Sol const& sol, double const& alpha, 
    int const& dis_flag, int const& times = 10000){
    vector<WVector> WVS = gene_all_rand_wvectors(D[0].size(), times, dis_flag);
    HP alpha_hit_times = 0;
    for(WVector const& wv: WVS){
        Score MaxS = 0;
        for(Tuple const& t: D) MaxS = max(MaxS, get_score(t, wv));
        for(TID const& tid: sol){
            if(get_score(D[tid], wv) >= MaxS * alpha){
                alpha_hit_times += 1;
                break;
            }
        }
    }
    return 1.0 * alpha_hit_times/times;
}