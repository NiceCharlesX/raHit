#include "MDRAH/io.cpp"
#include "MDRAH/basic.cpp"
#include "MDRAH/wvector.cpp"
#include "MDRAH/maxhp.cpp"
#include "MDRAH/test.cpp"

/* run ./main data n d r alpha dis_flag epsilon gamma
 * data : data name
 * n : number of tuple
 * d : dimension
 * r : output size
 * alpha : happiness ratio threshold of raHit
 * dis_flag : distribution flag
 * epsilon : absolute error of hit pro
 * gamma : discretization parameter
 */

int main(int argc, char *argv[])
{
    //input values
    string fileS = "anti4d";
    int n = 10000, d = 4, r = 5, dis_flag = 1, gamma = 12; //6 8 10 12 16 20
    double alpha = 0.99, epsilon = 0.03;
    
    //fixed values
    double delta = 0.01;

    ofstream RECORD("../record/mdrahp.csv", ios::app); //app = append
    if(argc > 1) fileS = string(argv[1]);
    if(argc > 2) n = stoi(string(argv[2]));
    if(argc > 3) d = stoi(string(argv[3]));
    if(argc > 4) r = stoi(string(argv[4]));
    if(argc > 5) alpha = stod(string(argv[5]));
    if(argc > 6) dis_flag = stoi(string(argv[6]));
    if(argc > 7) epsilon = stod(string(argv[7]));
    if(argc > 8) gamma = stoi(string(argv[8]));
    cout << "MDRAHP"<<","<<fileS<<",n="<<n<<",d="<<d<<",r="<<r<<",alpha="<<alpha<<",dflag=DF"<<dis_flag<<",epsilon="<<epsilon<<",gamma="<<gamma<<"\n";
    RECORD<<"MDRAHP"<<","<<fileS<<",n"<<n<<",d"<<d<<",r"<<r<<","<<alpha<<",DF"<<dis_flag<<","<<epsilon<<",g"<<gamma<<",";

    DataSet D;
    string fileName = "../data/"+ fileS +".data";
    read_file(D, n, d, fileName);

    struct timeval start, end;
    gettimeofday(&start,NULL);

    sort(D.begin(), D.end(), compare_tuples_by_sum());
    DataSet SD = compute_skyline(D);
    int sn = SD.size();
    
    int lambda = 1; //C(gamma+d-1, d) - C(gamma, d)
    for(int i = gamma + d - 1; i > gamma - 1; i--) lambda *= i;
    if(gamma >= d){
        int minus = 1;
        for(int i = gamma; i > gamma - d; i--) minus *= i;
        lambda -= minus;
    }
    for(int i = d; i > 0; i--) lambda /= i;

    int ssh = 0;
    if(lambda * log(2) > r * log(sn)){
        cout << "rlns\n";
        RECORD << "rlns,";
        ssh = (int)(log(2)+r*log(sn)+log(1.0/delta))/(2*pow(epsilon,2));
    }else{
        cout << "lambdaln2\n";
        RECORD << "lamln2,";
        ssh = (int)(log(2)+lambda*log(2)+log(1.0/delta))/(2*pow(epsilon,2));
    }

    map<WVector, int> WVC = gene_all_discretized_wvectors(d,ssh,dis_flag,gamma);
    map<ASet,int> ASC = compute_all_asets(SD,WVC,alpha);
    cout << "SD.size:" << sn << ", lambda:" << lambda << ", WVC.size:" << WVC.size() << ", ASC.size:" << ASC.size() << "\n";
    cout << "ssh:" << ssh << "\n";
    RECORD << "s" << sn << ",lam" << lambda << ","<< WVC.size() << "," << ASC.size() << ",m" << ssh << ",";

    int an = ASC.size();
    vector<ASet> ASV(an);
    vector<int> ACT(an); //from aid to corresponding count
    vector<set<AID>> TC(sn); //from sid to tuple's covering weight vectors
    map<SID, int> SCT; //from sid to the number of tuple's covering weight vectors
    initial_data_structures(ASC,ASV,ACT,TC,SCT);
    vector<bool> UC(an,true); //mark uncovered asets
    Sol sol = perform_greedy(r,ASV,ACT,TC,SCT,UC);
    
    gettimeofday(&end,NULL);

    HP pred_ahp = 0;
    for(AID aid = 0; aid < an; aid++) if(!UC[aid]) pred_ahp += ACT[aid];
    pred_ahp = pred_ahp/ssh;
    cout << "predicted alpha-hit pro:" << pred_ahp << "\n";
    RECORD << pred_ahp << ",";
    
    unsigned long timer = (end.tv_sec-start.tv_sec)*CONVERTER+(end.tv_usec-start.tv_usec);
    cout << fixed << setprecision(6) << "time: " << 1.0*timer/CONVERTER << "s\n";
    RECORD << fixed << setprecision(6) << 1.0*timer/CONVERTER << "s,";
    HP veri_ahp = test_solahp(SD, sol, alpha, dis_flag);
    cout << "verified alpha-hit pro:" << veri_ahp << "\n";
    RECORD << veri_ahp << endl;

    RECORD.close();
    cout << "\n";
    return 0;
}