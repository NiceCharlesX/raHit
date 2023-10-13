#include "struct.hpp"

/*
 * Generate a uniform random variable on top surface of [0,1]^d
 * i.e., assuming L_infinite = 1
 */
inline
WVector gene_a_wvector_in_linfinite(int const& d){
    WVector wv(d);
    int fixed_dim = lrand48()%d; //fix an attribute value to 1
    for(int i = 0; i < d; i++){
        if(i == fixed_dim) wv[i] = 1;
        else wv[i] = drand48();
    }
    return wv;
}

vector<WVector> gene_all_wvectors_in_linfinite(int const& d, int const& num){
    //for uniform distribution on [0,1]
    srand48((unsigned)time(nullptr)); 
    vector<WVector> WVS(num);
    for(int i = 0; i < num; i++) WVS[i] = gene_a_wvector_in_linfinite(d);
    return WVS;
}

/*
 * Generate a uniform random variable on \Delta^(d-1)
 * i.e., assuming L1 = 1
 */
inline
WVector gene_a_wvector_in_lone(int const& d){
    WVector wv(d);
    WVector mwv(d+1);
    mwv[0] = 0;
    for(int i = 1; i < d; i++) mwv[i] = drand48();
    mwv[d] = 1.0;
    sort(mwv.begin(), mwv.end());
    for(int i = 0; i < d; i++) wv[i] = mwv[i+1] - mwv[i];
    return wv;
}

vector<WVector> gene_all_wvectors_in_lone(int const& d, int const& num){
    //for uniform distribution on [0,1]
    srand48((unsigned)time(nullptr)); 
    vector<WVector> WVS(num);
    for(int i = 0; i < num; i++) WVS[i] = gene_a_wvector_in_lone(d);
    return WVS;
}

/*
 * Generate a uniform random variable on \DS^(d-1)
 * i.e., assuming L2 = 1
 */
inline
WVector gene_a_wvector_with_normal_dis(int const& d, mt19937_64 & gen, normal_distribution<double> & nd){
    WVector wv(d);
    for(int i = 0; i < d; i++) wv[i] = fabs(nd(gen));
    return wv;
}

vector<WVector> gene_all_wvectors_in_ltwo(int const& d, int const& num){
    //for normal distribution
    random_device rd; //Get a random seed from the OS entropy device, or whatever
    mt19937_64 gen(rd()); //Use the 64-bit Mersenne Twister 19937 generator
    normal_distribution<double> nd; //Standard normal distribution
    vector<WVector> WVS(num);
    for(int i = 0; i < num; i++) WVS[i] = gene_a_wvector_with_normal_dis(d,gen,nd);
    return WVS;
}

vector<WVector> gene_all_wvectors_with_normal_dis(int const& d, int const& num){
    //for normal distribution
    random_device rd; //Get a random seed from the OS entropy device, or whatever
    mt19937_64 gen(rd()); //Use the 64-bit Mersenne Twister 19937 generator
    normal_distribution<double> nd{0.5,0.25}; 
    vector<WVector> WVS(num);
    for(int i = 0; i < num; i++) WVS[i] = gene_a_wvector_with_normal_dis(d,gen,nd);
    return WVS;
}

//num is the number of random weightvectors
vector<WVector> gene_all_rand_wvectors(int const& d, int const& num, int const& dis_flag){
    switch (dis_flag){
    case 0: return gene_all_wvectors_in_linfinite(d, num);
    case 1: return gene_all_wvectors_in_lone(d, num); //candidate
    case 2: return gene_all_wvectors_in_ltwo(d, num);
    case 3: return gene_all_wvectors_with_normal_dis(d, num); //candidate
    default:
        cout<< "\nInput Parameter Error\n";
		exit(0);
    }
}
