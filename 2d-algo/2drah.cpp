#include "2DRAH/io.cpp"
#include "2DRAH/basic.cpp"
#include "2DRAH/dp.cpp"
#include "2DRAH/test.cpp"

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
    
    ofstream RECORD("../record/2drah.csv", ios::app); //app = append
    if(argc > 1) fileS = string(argv[1]);
    if(argc > 2) n = stoi(string(argv[2]));
    if(argc > 3) r = stoi(string(argv[3]));
    if(argc > 4) alpha = stod(string(argv[4]));
    cout << "2DRAH"<<","<<fileS<<",n="<<n<<",d=2,r="<<r<<",alpha="<<alpha<<"\n";
    RECORD<<"2DRAH"<<","<<fileS<<",n"<<n<<",d2,r"<<r<<","<<alpha<<",";

    DataSet D,SD; //SD: Sky(D)
    Sol sol; 
    SkyConvMap SCM;
    ConvSkyMap CSM;
    string fileName = "../data/"+ fileS +".data";
    read_file(D, n, fileName);

    struct timeval start, end;
    gettimeofday(&start,NULL);

    compute_sky_and_conv(D, SD, SCM, CSM);
    cout << "SD.size:"<<SD.size() << ", CSM.size:" << CSM.size() << "\n";
    RECORD << "s" << SD.size() << ",c" << CSM.size() << ",";
    if(alpha < 1){
        IntervalList IL = compute_IntervalList(SD, SCM, CSM, alpha); //compute_IntervalList_test(SD, SCM, CSM, alpha);
        if(IL.size() > r) sol = execute_dynamic_programming(IL, r);
        else{
            sol.hp = 1;
            for(Interval const& itv: IL) sol.sids.push_back(itv.sid);
        }
    }else{//alpha == 1
        sol = execute_greedy_process(SD,CSM,r);
    }

    for(SID const& sid: sol.sids) cout << "(" << SD[sid][0] << "," << SD[sid][1] << "),";
    cout << "\n";

    gettimeofday(&end,NULL);

    test_solhp(SD, sol, alpha);
    unsigned long timer = (end.tv_sec-start.tv_sec)*CONVERTER+(end.tv_usec-start.tv_usec);
    cout << fixed << setprecision(6) << "time: " << 1.0*timer/CONVERTER << "s\n";
    RECORD << fixed << setprecision(6) << 1.0*timer/CONVERTER << "s,";

    cout << "alpha-hit pro:" << sol.hp << "\n";
    RECORD << sol.hp << endl;
    RECORD.close();
    cout << "\n";
    return 0;
}