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
