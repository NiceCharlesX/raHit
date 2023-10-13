#include "struct.hpp"

struct compare_tuples_by_ycoord{ 
    bool operator() (Tuple const& a, Tuple const& b) const{
        if(a[1] == b[1]) return a[0] > b[0];
        return a[1] > b[1];
        
    }
};

struct  compare_intervals_by_left_border{
    bool operator() (Interval const& a, Interval const& b) const{
        if(a.p == b.p) return a.q > b.q; 
        return a.p < b.p;
        
    }
};

inline 
double cal_slope(Tuple const& a, Tuple const& b){
    return (a[1] - b[1])/(b[0] - a[0]); 
}

inline 
double cal_re_slope(Tuple const& a, Tuple const& b){ //reciprocal
    return (b[0] - a[0])/(a[1] - b[1]); 
}

//compute skyline SD and label tuples in ConvD
void compute_sky_and_conv(DataSet & D, DataSet & SD, SkyConvMap & SCM, ConvSkyMap & CSM){
    sort(D.begin(), D.end(), compare_tuples_by_ycoord());//sort tuples by the y-coordinate in the descending order
    Coord MaxX = -1;//largest x-coordinate
    for(Tuple const & t: D){
        if(t[0] > MaxX){ //is in skyline
            SD.push_back(t);
            MaxX = t[0];
        }
    }
    //cout << "SkyD.size:" << SD.size() << ", ";
    
    CSM.push_back(0); 
    if(SD.size() > 1){ 
        CSM.push_back(1); 
        CID cid = 1; 
        if(SD.size() > 2){
            for(SID sid = 2; sid < SD.size(); sid++){
                while(cid > 0 && 
                    cal_slope(SD[CSM[cid - 1]], SD[CSM[cid]]) >= cal_slope(SD[CSM[cid]], SD[sid])){
                    
                    CSM.pop_back();
                    cid --;
                }
                CSM.push_back(sid); 
                cid ++;
            }
        }
    }
    //cout << "ConvD.size:" << CSM.size() << "\n";

    SCM = vector<CID>(SD.size(), -1);
    for(CID cid = 0; cid < CSM.size(); cid++){
        SCM[CSM[cid]] = cid;
    }
}

//Assumption: alpha < 1
IntervalList compute_IntervalList(DataSet const& SD, SkyConvMap const& SCM, ConvSkyMap const& CSM, 
    double const& alpha){
    IntervalList IL;

    Coord MaxX = SD[SD.size()-1][0], MaxY = SD[0][1];
    CID cid = 0; //cid of latest vertex/tuple in ConvD 
    for(SID sid = 0; sid < SD.size(); sid++){
        if(SCM[sid] >= 0) cid = SCM[sid]; 
        Tuple t = {SD[sid][0]/alpha, SD[sid][1]/alpha}; 

        if(cid + 1 < CSM.size() && t[1] - SD[CSM[cid+1]][1] <= 
            (SD[CSM[cid+1]][0] - t[0]) * cal_slope(SD[CSM[cid]], SD[CSM[cid+1]])){
            //t is not in Conv(D+t)
            
            continue;
        }

        double p, q; //left/right border of interval
        if(t[1] >= MaxY) p = 0;
        else{
            CID start = 0, end = cid;
            while(end > start){
                CID mid = start + (end - start)/2;
                double slope_diff = cal_slope(SD[CSM[mid]], t) - cal_slope(SD[CSM[mid + 1]], t);
                
                if(slope_diff > 0) end = mid;
                else if(slope_diff < 0) start = mid + 1;
                else{
                    start = mid; 
                    end = mid;
                }
            }
            p = (SD[CSM[start]][1] - t[1]) / (SD[CSM[start]][1] - t[1] + t[0] - SD[CSM[start]][0]); 
        }

        if(t[0] >= MaxX) q = 1;
        else{
            CID start = cid + 1, end = CSM.size() - 1;
            while(end > start){
                CID mid = start + (end - start)/2;
                double re_slope_diff = cal_re_slope(t, SD[CSM[mid]]) - cal_re_slope(t, SD[CSM[mid + 1]]);
                
                if(re_slope_diff > 0) end = mid;
                else if(re_slope_diff < 0) start = mid + 1;
                else{
                    start = mid + 1;
                    end = mid + 1;
                }
            }
            q = (t[1] - SD[CSM[start]][1]) / (t[1] - SD[CSM[start]][1] + SD[CSM[start]][0] - t[0]);
        }
        IL.push_back({sid, p, q});
        //cout << "(" << SD[sid][0] << "," << SD[sid][1] << "), [" << p << "," << q << "]\n";
    }

    cout << "IL.size:" << IL.size() << ", ";
    sort(IL.begin(), IL.end(), compare_intervals_by_left_border());//sort intervals by left border in increasing order

    IntervalList Final_IL;
    Final_IL.push_back(IL[0]);
    if(IL.size() > 1){
        for(int i = 1; i < IL.size(); i++){
            //cout << "[" << IL[i].p << "," << IL[i].q << "]\n";
            if(IL[i].q > Final_IL.back().q) Final_IL.push_back(IL[i]);
        }
    }
    
    cout << "FIL.size:" << Final_IL.size() << "\n";
    return Final_IL;
}


//Assumption: alpha < 1
IntervalList compute_IntervalList_test(DataSet const& SD, SkyConvMap const& SCM, ConvSkyMap const& CSM, 
    double const& alpha){
    IntervalList IL;

    Coord MaxX = SD[SD.size()-1][0], MaxY = SD[0][1];
    CID cid = 0; //cid of latest vertex/tuple in ConvD 
    for(SID sid = 0; sid < SD.size(); sid++){
        if(SCM[sid] >= 0) cid = SCM[sid]; 
        Tuple t = {SD[sid][0]/alpha, SD[sid][1]/alpha}; 

        if(cid + 1 < CSM.size() && t[1] - SD[CSM[cid+1]][1] <= 
            (SD[CSM[cid+1]][0] - t[0]) * cal_slope(SD[CSM[cid]], SD[CSM[cid+1]])){
            //t is not in Conv(D+t)
            
            continue;
        }

        double p, q; //left/right border of interval
        if(t[1] >= MaxY) p = 0;
        else{
            CID start = cid;
            double old_slope = cal_slope(SD[CSM[start]], t), new_slope;
            while(start > 0){
                new_slope = cal_slope(SD[CSM[start-1]], t);
                if(new_slope < old_slope) break;
                old_slope = new_slope;
                start --;
            }
            p = (SD[CSM[start]][1] - t[1]) / (SD[CSM[start]][1] - t[1] + t[0] - SD[CSM[start]][0]); 
        }

        if(t[0] >= MaxX) q = 1;
        else{
            CID start = cid + 1;
            double old_re_slope = cal_re_slope(t, SD[CSM[start]]), new_re_slope;
            while(start < CSM.size()-1){
                new_re_slope = cal_re_slope(t, SD[CSM[start+1]]);
                if(new_re_slope < old_re_slope) break;
                old_re_slope = new_re_slope;
                start ++;
            }
            q = (t[1] - SD[CSM[start]][1]) / (t[1] - SD[CSM[start]][1] + SD[CSM[start]][0] - t[0]);
        }
        IL.push_back({sid, p, q});
        //cout << "(" << SD[sid][0] << "," << SD[sid][1] << "), [" << p << "," << q << "]\n";
    }

    cout << "IL.size:" << IL.size() << ", ";
    sort(IL.begin(), IL.end(), compare_intervals_by_left_border());//sort intervals by left border in increasing order

    IntervalList Final_IL;
    Final_IL.push_back(IL[0]);
    if(IL.size() > 1){
        for(int i = 1; i < IL.size(); i++){
            //cout << "[" << IL[i].p << "," << IL[i].q << "]\n";
            if(IL[i].q > Final_IL.back().q) Final_IL.push_back(IL[i]);
        }
    }
    
    cout << "FIL.size:" << Final_IL.size() << "\n";
    return Final_IL;
}


struct  compare_intervals_by_extent{ //for priority_queue
    bool operator() (Interval const& a, Interval const& b) const{
        return a.q - a.p > b.q - b.p; 
    }
    
};

//Assumption: alpha = 1
Sol execute_greedy_process(DataSet const& SD, ConvSkyMap const& CSM, int const& r){
    priority_queue<Interval, vector<Interval>, compare_intervals_by_extent> MaxRI;
    for(CID cid = 0; cid < CSM.size(); cid++){
        double p, q; //left/right border of interval
        if(cid > 0) p = (SD[CSM[cid-1]][1] - SD[CSM[cid]][1]) / (SD[CSM[cid-1]][1] - SD[CSM[cid]][1] + SD[CSM[cid]][0] - SD[CSM[cid-1]][0]); 
        else p = 0;

        if(cid < CSM.size() - 1) q = (SD[CSM[cid]][1] - SD[CSM[cid+1]][1]) / (SD[CSM[cid]][1] - SD[CSM[cid+1]][1] + SD[CSM[cid+1]][0] - SD[CSM[cid]][0]);
        else q = 1;

        if(MaxRI.size() < r) MaxRI.push({CSM[cid], p, q});
        else if(MaxRI.top().q - MaxRI.top().p < q - p){
            MaxRI.pop();
            MaxRI.push({CSM[cid], p, q});
        }
        //cout << "(" << SD[CSM[cid]][0] << "," << SD[CSM[cid]][1] << "), [" << p << "," << q << "], extent:" << q - p << "\n";
    }

    Sol sol;
    sol.hp = 0;
    while(!MaxRI.empty()){
        Interval itv = MaxRI.top();
        sol.hp += itv.q - itv.p;
        sol.sids.push_back(itv.sid);
        MaxRI.pop();
    }
    return sol;
}