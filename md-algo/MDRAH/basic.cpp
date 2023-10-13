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
        Coord sa = 0, sb = 0; //sums
        for(int i = 0; i < a.size(); i++){
            sa += a[i];
            sb += b[i];
        }
        return sa > sb;
    }
};

/*
 * determine whether a is dominated by b
 */
bool is_dominated(Tuple const& a, Tuple const & b){
    for(int i = 0; i < a.size(); i++){
        if(a[i] > b[i]) return false;
    }
    return true;
}

DataSet compute_skyline(DataSet const& D){
    DataSet SD;
    for(Tuple const& t: D){
        bool notdom = true;
        for(Tuple st: SD){
            if(is_dominated(t, st)){
                notdom = false;
                break;
            }
        }
        if(notdom) SD.push_back(t);
    }
    return SD;
}
