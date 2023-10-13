#include "MDRKH/io.cpp"
#include "MDRKH/basic.cpp"
#include "MDRKH/wvector.cpp"
#include "MDRKH/skyline.cpp"
#include "MDRKH/maxuhp.cpp"
#include "MDRKH/test.cpp"

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
    double delta = 0.01, epsilon = 0.03;
    int k = 16; //4, 6, 8, 10, 12, 16

    ofstream RECORD("../record/mdrkh.csv", ios::app); //app = append
    if(argc > 1) fileS = string(argv[1]);
    if(argc > 2) n = stoi(string(argv[2]));
    if(argc > 3) d = stoi(string(argv[3]));
    if(argc > 4) r = stoi(string(argv[4]));
    if(argc > 5) alpha = stod(string(argv[5]));
    if(argc > 6) dis_flag = stoi(string(argv[6]));
    cout << "MDRKH"<<","<<fileS<<",n="<<n<<",d="<<d<<",r="<<r<<",alpha="<<alpha<<",dflag=DF"<<dis_flag<<",k="<<k<<"\n";
    RECORD<<"MDRKH"<<","<<fileS<<",n"<<n<<",d"<<d<<",r"<<r<<","<<alpha<<",DF"<<dis_flag<<",k"<<k<<",";

    DataSet D;
    string fileName = "../data/"+ fileS +".data";
    read_file(D, n, d, fileName); //read_4file(D, n, d, fileName);

    struct timeval start, end;
    gettimeofday(&start,NULL);

    sort(D.begin(), D.end(), compare_tuples_by_sum());
    vector<DataSet> LD = compute_skyband(D,k);
    DataSet SD = LD[0],kSD; //SD: Sky(D); kSD: k-Sky(D)/Sky(D)
    for(int i = 1; i < k; i++) kSD.insert(kSD.end(), LD[i].begin(), LD[i].end());
    int sn = SD.size(), ssh = (int)(log(2)+r*log(sn)+log(1.0/delta))/(2*pow(epsilon,2));

    vector<WVector> WVS = gene_all_wvectors(d,ssh,dis_flag);
    map<KSet,KI> KSS = compute_all_ksets(SD,kSD,WVS,k);
    cout << "SD.size:" << sn << ", kSD.size:" << kSD.size() << ", KSS.size:" << KSS.size() << "\n";
    cout << "ssh:" << ssh << "\n";
    RECORD << "s" << sn << ",ks" << kSD.size() << "," << KSS.size() << ",m" << ssh << ",";

    int kn = KSS.size();
    vector<KSet> KSV(kn);
    vector<WVector> CWV(kn); //from kid to corresponding weight vector
    vector<int> KCT(kn); //from kid to corresponding count
    vector<set<KID>> TC(sn); //from sid to tuple's covering weight vectors
    map<SID, int> SCT; //from sid to the number of tuple's covering weight vectors
    initial_data_structures(KSS,KSV,CWV,KCT,TC,SCT);
    vector<bool> UC(kn,true); //mark uncovered ksets
    Sol sol = perform_greedy(r,KSV,KCT,TC,SCT,UC);
    
    int round = 0;
    Sol fsol = perform_clustering(SD,r,sol,round,CWV,KCT,TC);
    cout << "round num:" << round << ", ";
    RECORD << round << ",";

    gettimeofday(&end,NULL);

    HP k_hp = test_solkhp(SD, kSD, fsol, k, dis_flag);
    cout << "k-hit pro:" << k_hp << "\n";
    RECORD << "khp" << k_hp << ",";

    unsigned long timer = (end.tv_sec-start.tv_sec)*CONVERTER+(end.tv_usec-start.tv_usec);
    cout << fixed << setprecision(6) << "time: " << 1.0*timer/CONVERTER << "s\n";
    RECORD << fixed << setprecision(6) << 1.0*timer/CONVERTER << "s,";
    HP alpha_hp = test_solahp(SD, fsol, alpha, dis_flag);
    cout << "alpha-hit pro:" << alpha_hp << "\n";
    RECORD << alpha_hp << endl; 

    RECORD.close();
    cout << "\n";
    return 0;
}