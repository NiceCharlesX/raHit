#include "HDRRM/io.cpp"
#include "HDRRM/basic.cpp"
#include "HDRRM/wvector.cpp"
#include "HDRRM/sort.cpp"
#include "HDRRM/sc.cpp"
#include "HDRRM/test.cpp"

/* run ./main data n d r alpha dis_flag
 * data : data name
 * n : number of tuple
 * d : dimension
 * r : output size
 * alpha : happiness ratio threshold of raHit
 * dis_flag : distribution flag
 */

int main(int argc, char *argv[]){
    //input values
    string fileS = "anti4d";
    int n = 10000, d = 4, r = 5, dis_flag = 1;
    double alpha = 0.99;

    //fixed values
    double delta = 0.03; //absolute error
    int gamma = 6;
    
    ofstream RECORD("../record/hdrrm.csv", ios::app); //app = append
    if(argc > 1) fileS = string(argv[1]);
    if(argc > 2) n = stoi(string(argv[2]));
    if(argc > 3) d = stoi(string(argv[3]));
    if(argc > 4) r = stoi(string(argv[4]));
    if(argc > 5) alpha = stod(string(argv[5]));
    if(argc > 6) dis_flag = stoi(string(argv[6]));
    cout << "HDRRM"<<","<<fileS<<",n="<<n<<",d="<<d<<",r="<<r<<",alpha="<<alpha<<",dflag=DF"<<dis_flag<<",gamma="<<gamma<<"\n";
    RECORD<<"HDRRM"<<","<<fileS<<",n"<<n<<",d"<<d<<",r"<<r<<","<<alpha<<",DF"<<dis_flag<<",g"<<gamma<<",";

    DataSet D;
    string fileName = "../data/"+ fileS +".data";
    read_file(D, n, d, fileName);

    struct timeval start, end;
    gettimeofday(&start,NULL);

    sort(D.begin(), D.end(), compare_tuples_by_sum());
    Sol basis = get_basis(D,d), sol;
    r = r + basis.size(); //The size bound of the solution does not consider basis

    filter_dataset(d, D, n, basis);
    int ssh = (int)((r-d)*log(n-d)+log(n-r+1)+log(n))/(2*pow(delta-1.0/n,2));
    int sd = pow(gamma+1,d-1);
    cout << "datasize after filtering:" << n << "\n";
    cout << "ssh:" << ssh << ", sd:" << sd << "\n";
    RECORD << n << ",m" << ssh << "," << sd << ",";

    Rank k = 1;
    vector<vector<SItem>> SILS(ssh+sd);
    vector<WVector> WVS = gene_all_rand_wvectors(d, ssh, dis_flag);
    add_non_rand_wvectors(WVS,d,gamma);
    sol = set_cover(init_sitem_lists(SILS, WVS, D, basis), n);

    while(sol.size() > r - basis.size()){
        k *= 2;
        sol = set_cover(incre_compute_ksets(SILS,k), n);
    }//double
    Rank low = k/2+1, high = k, fk = k;
    while(low < high){
        k = (low + high)/2;
        Sol msol = set_cover(decre_compute_ksets(SILS, k, low, high), n); 
        if(msol.size() + basis.size() <= r){
            high = k;
            sol = msol;
            fk = k;
        }else low = k + 1;
    }//binary search
    sol.insert(basis.begin(), basis.end());
    cout << "rank-regret:" << fk << ", sol.size:" << sol.size() << "\n";
    RECORD << "rr" << fk << "," << sol.size() << ",";

    gettimeofday(&end,NULL);

    unsigned long timer = (end.tv_sec-start.tv_sec)*CONVERTER+(end.tv_usec-start.tv_usec);
    cout << fixed << setprecision(6) << "total time: " << 1.0*timer/CONVERTER << "s\n";
    RECORD << fixed << setprecision(6) << 1.0*timer/CONVERTER << "s,";
    HP alpha_hp = test_solahp(D, sol, alpha, dis_flag);
    cout << "alpha-hit pro:" << alpha_hp << "\n";
    RECORD << alpha_hp << endl;

    RECORD.close();
    cout << "\n";
    return 0;
}