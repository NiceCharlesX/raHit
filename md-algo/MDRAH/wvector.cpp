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
        if(i == fixed_dim) wv[i] = 1.0;
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
    mwv[0] = 0.0;
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
vector<WVector> gene_all_wvectors(int const& d, int const& num, int const& dis_flag){
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

inline
WVector transform_randwv_to_discwv_linfinite(WVector const& rand_wv,
    int const& gamma){
    int max_i = max_element(rand_wv.begin(),rand_wv.end()) - rand_wv.begin(); 
/*
    int max_dim = 0; 
    Weight max_wei = rand_wv[0];
    for(int i = 1; i < rand_wv.size(); i++){
        if(rand_wv[i] > max_wei){
            max_dim = i;
            max_wei = rand_wv[i];
        }
    }
*/
    WVector disc_wv(rand_wv.size(), 0.0);
    for(int i = 0; i < rand_wv.size(); i++){
        if(i == max_i) disc_wv[i] = 1.0;
        else{
            disc_wv[i] = floor(rand_wv[i]/rand_wv[max_i] * gamma) / gamma;
            if(disc_wv[i] < 1.0) disc_wv[i] += 0.5/gamma;
            else{ 
                cout << "Extreme Case!!!\n";
                disc_wv[i] -= 0.5/gamma;
            }
        }
    }
    return disc_wv;
}

inline
WVector transform_randwv_to_discwv_lone_backup(WVector const& rand_wv,
    int const& gamma){
    Weight sum_randw = accumulate(rand_wv.begin(), rand_wv.end(), 0.0); 
    WVector disc_wv(rand_wv.size(), 0.0);
    Weight sum_discw = 0.0;

    for(int i = 0; i < rand_wv.size()-1; i++){ 
        disc_wv[i] = floor(rand_wv[i]/sum_randw * gamma) / gamma;
        sum_discw += disc_wv[i];
    }
    disc_wv[rand_wv.size()-1] = 1.0 - sum_discw;

    return disc_wv;
}

inline
WVector transform_randwv_to_discwv_lone(WVector const& rand_wv,
    int const& gamma){
    Weight sum_randw = accumulate(rand_wv.begin(), rand_wv.end(), 0.0); 
    WVector disc_wv(rand_wv.size(), 0.0);
    Weight sum_discw = 0.0;

    for(int i = 0; i < rand_wv.size(); i++){ 
        disc_wv[i] = floor(rand_wv[i]/sum_randw * gamma) / gamma;
        sum_discw += disc_wv[i];
    }
    if(sum_discw < 1.0){
        for(int i = 0; i < rand_wv.size(); i++){
            disc_wv[i] += (1.0 - sum_discw)/rand_wv.size();
        }
    }else{ 
        cout << "Extreme Case!!!\n";
        for(int i = 0; i < rand_wv.size(); i++){
            if(disc_wv[i] > 0){
                disc_wv[i] -= 1.0/gamma;
                break;
            }
        }
        for(int i = 0; i < rand_wv.size(); i++){
            disc_wv[i] += 1.0/gamma/rand_wv.size();
        }
    }

    return disc_wv;
}

map<WVector, int> gene_all_discretized_wvectors(int const& d, int const& num, int const& dis_flag,
    int const& gamma){
    map<WVector, int> WVC;
    vector<WVector> WVS = gene_all_wvectors(d, num, dis_flag);
    for(WVector const rwv: WVS){
        WVector dwv = transform_randwv_to_discwv_lone(rwv, gamma);
        //WVector dwv = transform_randwv_to_discwv_lone_backup(rwv, gamma);
        //WVector dwv = transform_randwv_to_discwv_linfinite(rwv, gamma);
        if(WVC.count(dwv) == 0)
            WVC.insert({dwv, 1});
        else WVC[dwv] += 1;
    } 

    /*for(auto const& mpair: WVC){
        cout << "(";
        for(int i = 0; i < d; i++) cout << mpair.first[i] << ",";
        cout << mpair.second << ")\n";
    }*/
    return WVC;
}