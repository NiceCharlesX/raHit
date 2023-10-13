#include "MDRMS/io.cpp"
#include "MDRMS/wvector.cpp"
#include "MDRMS/basic.cpp"
#include "MDRMS/sc.cpp"
#include "MDRMS/test.cpp"

/* run ./main data n d r alpha dis_flag
 * data : data name
 * n : number of tuple
 * d : dimension
 * r : output size
 * alpha : happiness ratio threshold of raHit
 * dis_flag : distribution flag
 */

int main(int argc, char *argv[])
{
    //input values
    string fileS = "anti4d";
    int n = 10000, d = 4, r = 5, dis_flag = 1;
    double alpha = 0.99;
    
    //fixed values
    int gamma = 6;
    
    ofstream RECORD("../record/mdrms.csv", ios::app); //app = append
    if(argc > 1) fileS = string(argv[1]);
    if(argc > 2) n = stoi(string(argv[2]));
    if(argc > 3) d = stoi(string(argv[3]));
    if(argc > 4) r = stoi(string(argv[4]));
    if(argc > 5) alpha = stod(string(argv[5]));
    if(argc > 6) dis_flag = stoi(string(argv[6]));
    cout << "MDRMS"<<","<<fileS<<",n="<<n<<",d="<<d<<",r="<<r<<",alpha="<<alpha<<",dflag=DF"<<dis_flag<<",gamma="<<gamma<<"\n";
    RECORD<<"MDRMS"<<","<<fileS<<",n"<<n<<",d"<<d<<",r"<<r<<","<<alpha<<",DF"<<dis_flag<<",g"<<gamma<<",";

    DataSet D;
    string fileName = "../data/"+ fileS +".data";
    read_file(D, n, d, fileName);

    struct timeval start, end;
    gettimeofday(&start,NULL);

    sort(D.begin(), D.end(), compare_tuples_by_sum());
    DataSet SD = get_skyline(D); //compute skyline
    int sn = SD.size(), sd = pow(gamma+1,d-1);

    vector<vector<RRatio>> MAT(sn, vector<RRatio>(sd));
    vector<RRatio> RRL;
    vector<WVector> WVS = gene_all_non_rand_wvectors(d, gamma);
    init_rratio_matrix(MAT,RRL,WVS,SD);
    cout << "SD.size:" << sn << ", WVS.size:" << WVS.size() << "\n";
    RECORD << "s" << sn << "," << WVS.size() << ","; 

    Sol sol;
    int low = 0, high = RRL.size() - 1;
    RRatio frr;
    while(low < high){
        int mid = (low + high)/2;
        Sol msol = mrst_oracle(MAT, RRL[mid], sn, sd);
        if(msol.size() <= r){
            sol = msol;
            frr = RRL[mid];
            high = mid - 1;
        }else low = mid + 1;
    }//binary search
    cout << "maximum regret-ratio:" << frr << "\n";
    RECORD << frr << ",";

    gettimeofday(&end,NULL);

    unsigned long timer = (end.tv_sec-start.tv_sec)*CONVERTER+(end.tv_usec-start.tv_usec);
    cout << fixed << setprecision(6) << "total time: " << 1.0*timer/CONVERTER << "s\n";
    RECORD << fixed << setprecision(6) << 1.0*timer/CONVERTER << "s,";
    HP alpha_hp = test_solahp(SD, sol, alpha, dis_flag);
    cout << "alpha-hit pro:" << alpha_hp << "\n";
    RECORD << alpha_hp << endl;

    RECORD.close();
    cout << "\n";
    return 0;
}