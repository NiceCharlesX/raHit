#include "struct.hpp"

inline
Score get_score(Tuple const& t, WVector const& wv){
    Score re = 0;
    for(int i = 0; i < t.size(); i++) re += wv[i] * t[i];
    return re;
}

/*
 * Divide the d-1-dimensional polar coordinate space into gamma^(d-1) lattices
 * to obtain (gamma+1)^(d-1) boundary weight vectors.
 */
void add_non_rand_wvectors(vector<WVector> & WVS, int const& d, int const& gamma){
    int sd = pow(gamma+1, d-1);

    vector<int> ivector(d-1,0);  //Intermediate vector
    for(int j = 0; j < sd; j++){
        WVector wv(d);
        Weight ivalue = 1.0; //Intermediate value
        for(int i = d-1; i > 0; i--){
            wv[i] = ivalue * cos(ivector[i-1] * M_PI/2.0/gamma);
            ivalue = ivalue * sin(ivector[i-1] * M_PI/2.0/gamma);
        }
        wv[0] = ivalue;

        WVS.push_back(wv);

        for(int i = 0; i < d-1; i++){
            if(ivector[i]<gamma){
                ivector[i] += 1;
                break;
            }else ivector[i] = 0;
        }
    }
}

Sol get_basis(DataSet const& D, int const& d){
    vector<Coord> maxvs(d,0); //maximum value on each dimension
    vector<TID> ctids(d,-1); //corresponding point id on each dimension
    for(int i = 0; i < D.size(); i++){
        for(int j = 0; j < d; j++){
            if(D[i][j] > maxvs[j]){
                maxvs[j] = D[i][j];
                ctids[j] = i;
            }
        }
    }
    return Sol(ctids.begin(), ctids.end());
}

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
 * Remove tuples in D dominated by basis
 */
void filter_dataset(int const& d, 
    DataSet & D, int & n, Sol & basis){

    DataSet FD; //filtered dataset
    Sol fbasis;
    TID count = 0;
    for(TID const& tid: basis){
        FD.push_back(D[tid]);
        fbasis.insert(count);
        count++;
    }
    for(Tuple const& t:D){
        bool notdom = true;
        for(TID const& tid: fbasis){
            if(is_dominated(t, FD[tid])) notdom = false;
        }
        if(notdom) FD.push_back(t);
    }
    D = FD;
    n = D.size();
    basis = fbasis;
}
