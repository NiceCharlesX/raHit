#include "2DRKH/io.cpp"
#include "2DRKH/basic.cpp"
#include "2DRKH/test.cpp"

/* run ./main data n r alpha
 * data : data name
 * n : number of tuple
 * r : output size
 * alpha : threshold
 */

int main(int argc, char *argv[]){
    //input values
    string fileS = "anti2d";
    int n = 10000, r = 5;
    double alpha = 0.999;

    //fixed values
    int k = 16; 

    ofstream RECORD("../record/2drkh.csv", ios::app); //app = append
    if(argc > 1) fileS = string(argv[1]);
    if(argc > 2) n = stoi(string(argv[2]));
    if(argc > 3) r = stoi(string(argv[3]));
    if(argc > 4) alpha = stod(string(argv[4]));
    cout << "2DRKH"<<","<<fileS<<",n="<<n<<",d=2,r="<<r<<",alpha="<<alpha<<",k="<<k<<"\n";
    RECORD<<"2DRKH"<<","<<fileS<<",n"<<n<<",d2,r"<<r<<","<<alpha<<",k"<<k<<",";

    DataSet D,SD,kSD; //SD: Sky(D), kSD: k-Sky(D)/Sky(D)
    SkyRankMap SRM;
    InterSQ ISQ;
    vector<vector<Sol>> MAT;
    string fileName = "../data/"+ fileS +".data";
    read_file(D, n, fileName);

    struct timeval start, end;
    gettimeofday(&start,NULL);

    compute_skyband(D, k, SD, kSD, SRM); //at the same, initial the rank of each skyline tuple
    cout << "SD.size:" << SD.size() << ", kSD.size:" << kSD.size() << "\n";
    RECORD << "s" << SD.size()<<",ks"<<kSD.size()<<",";
    int s = SD.size(); //skyline size
    get_intersq(SD, kSD, ISQ);
    init_matrix(s,r,MAT);
    Sol sol = process_intersq(s, r, k, MAT, SRM, ISQ);
    for(SID const& sid: sol.sids) cout << "(" << SD[sid][0] << "," << SD[sid][1] << "),";
    cout << "\nk-hit pro:" << sol.hp << ", ";
    RECORD << "khp" << sol.hp << ",";

    gettimeofday(&end,NULL);

    test_skyrankmap(SRM);
    HP alpha_hp = test_solhp(SD, sol, alpha);
    unsigned long timer = (end.tv_sec-start.tv_sec)*CONVERTER+(end.tv_usec-start.tv_usec);
    cout << fixed << setprecision(6) << "time: " << 1.0*timer/CONVERTER << "s\n";
    RECORD << fixed << setprecision(6) << 1.0*timer/CONVERTER << "s,";
    cout << "alpha-hit pro:" << alpha_hp << "\n";
    RECORD<<alpha_hp<<endl;

    RECORD.close();
    cout << "\n";
    return 0;
}