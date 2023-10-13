#include "struct.hpp"

void test_solhp(DataSet const& SD, Sol const& sol, 
    double const& alpha, int const& times = 10000){
    srand48((unsigned)time(nullptr));
    int alpha_hit_times = 0, one_hit_times = 0;
    for(int i = 0; i < times; i++){
        Coord w0 = drand48();
        Score MaxS = 0, MaxS_sol = 0;
        for(SID const& sid: sol.sids) MaxS_sol = max(MaxS_sol, SD[sid][0] * w0 + SD[sid][1] * (1 - w0));
        for(Tuple const& t: SD) MaxS = max(MaxS, t[0] * w0 + t[1] * (1 - w0));

        if(MaxS_sol >= MaxS * alpha) alpha_hit_times++;
        if(MaxS_sol >= MaxS) one_hit_times++;
    }
    cout << "delta:" << sol.hp - 1.0 * alpha_hit_times/times << ", ";
    cout << "1-hit pro:" << 1.0 * one_hit_times/times << "\n";
}