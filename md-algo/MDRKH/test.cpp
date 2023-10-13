#include "struct.hpp"

HP test_solahp(DataSet const& SD, Sol const& sol, double const& alpha, 
    int const& dis_flag, int const& times = 10000){
    vector<WVector> WVS = gene_all_wvectors(SD[0].size(), times, dis_flag);
    HP alpha_hit_times = 0;
    for(WVector const& wv: WVS){
        Score MaxS = 0;
        for(Tuple const& t: SD) MaxS = max(MaxS, get_score(t, wv));
        for(SID const& sid: sol){
            if(get_score(SD[sid], wv) >= MaxS * alpha){
                alpha_hit_times += 1;
                break;
            }
        }
    }
    return 1.0 * alpha_hit_times/times;
}

HP test_solkhp(DataSet const& SD, DataSet const& kSD, Sol const& sol, int const& k, 
    int const& dis_flag, int const& times = 10000){
    vector<WVector> WVS = gene_all_wvectors(SD[0].size(), times, dis_flag);
    HP k_hit_times = 0;
    for(WVector const& wv: WVS){
        Score MaxS_Sol = 0;
        for(SID const& sid: sol) MaxS_Sol = max(MaxS_Sol, get_score(SD[sid], wv));
        Rank rank = 1;
        for(Tuple const& t: SD){
            if(get_score(t, wv) > MaxS_Sol) rank++;
        }
        for(Tuple const& t: kSD){
            if(get_score(t, wv) > MaxS_Sol) rank++;
        }
        if(rank <= k) k_hit_times++;
    }
    return 1.0 * k_hit_times/times;
}