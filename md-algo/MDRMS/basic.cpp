#include "struct.hpp"

inline
Score get_score(Tuple const& t, WVector const& wv){
    Score re = 0;
    for(int i = 0; i < t.size(); i++) re += wv[i] * t[i];
    return re;
}

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

/*
 * Divide the d-1-dimensional polar coordinate space into gamma^(d-1) lattices
 * to obtain (gamma+1)^(d-1) boundary weight vectors.
 */
vector<WVector> gene_all_non_rand_wvectors(int const& d, int const& gamma){ //utility vector num
    int sd = pow(gamma + 1, d - 1);
    vector<WVector> WVS(sd);
    vector<int> ivector(d-1,0); //intermediate vector
    for(int j = 0; j < sd; j++){
        WVector wv(d);
        Weight ivalue = 1.0; //intermediate value
        for(int i = d - 1; i > 0; i--){
            wv[i] = ivalue * cos(ivector[i-1] * M_PI/2.0/gamma);
            ivalue = ivalue * sin(ivector[i-1] * M_PI/2.0/gamma);
        }
        wv[0] = ivalue;
        WVS[j] = wv;
        for(int i = 0; i < d-1; i++){
            if(ivector[i]<gamma){
                ivector[i] += 1;
                break;
            }else ivector[i] = 0;
        }
    }
    return WVS;
}

/*
 * sn is the size of SD (Skyline)
 * MAT[i][j], i denots SID, j denotes WID 
 */
void init_rratio_matrix(vector<vector<RRatio>> & MAT, vector<RRatio> & RRL, 
    vector<WVector> const& WVS, DataSet const& SD){
    for(WID j = 0; j < WVS.size(); j++){
        Score maxs = 0; //maximum score in SD
        for(Tuple t: SD) maxs = max(maxs, get_score(t, WVS[j]));
        for(SID i = 0; i < SD.size(); i++){
            RRatio rr = (maxs - get_score(SD[i], WVS[j])) / maxs;
            MAT[i][j] = rr;
            if(rr > 0) RRL.push_back(rr);
        }
    }
    sort(RRL.begin(), RRL.end()); //in the ascending order
}

//True: a is dominated by b;
bool is_dominated(Tuple const& a, Tuple const & b){
    for(int i = 0; i < a.size(); i++){
        if(a[i] > b[i]) return false;
    }
    return true;
}

DataSet get_skyline(DataSet const& D){
    DataSet SD;
    for(Tuple const& t: D){
        bool notdom = true;
        for(Tuple st: SD){
            if(is_dominated(t, st)) notdom = false;
        }
        if(notdom) SD.push_back(t);
    }
    return SD;
}