#include "struct.hpp"

HP get_penalty(Rank const& rank, int const& k){
    if(rank <= k) return 1;//1/rank; //1/log2(1+rank);
    else return 0;
}
HP get_cumul_pr(Coord const& start, Coord const& end){
    return end - start;
}

struct compare_tuples_by_ycoord{
    bool operator() (Tuple const& a, Tuple const& b) const{
        if(a[1] == b[1]) return a[0] > b[0];
        return a[1] > b[1];
    }
};

//compute skyline SD, compute k-skyband SD\cup kSD. 
void compute_skyband(DataSet & D, int const& k, 
    DataSet & SD, DataSet & kSD, SkyRankMap & SRM){

    sort(D.begin(), D.end(), compare_tuples_by_ycoord());//sort tuples by the y-coordinate in the descending order
    Coord MaxX = 0;//largest x-coordinate
    priority_queue<Coord, vector<Coord>, greater<Coord>> MaxKX;//k largest x-coordinate, min-heap
    Rank rank = 0;
    for(Tuple const & t: D){
        if(MaxKX.size() < k || t[0] > MaxKX.top()){// is in kskyband, if t[0] <= MaxKX.top(), there must be k tuples dominating t
            rank ++;
            if(t[0] > MaxX){ //is in skyline
                SD.push_back(t);
                SRM.push_back(rank); //initial rank of each skyline tuple
                MaxX = t[0];
            }else kSD.push_back(t);

            if(MaxKX.size() == k) MaxKX.pop(); //.top() is the k-th largest
            else if(MaxKX.size() > k){
                cout << "\nerror 6!\n";
                abort();
            }

            MaxKX.push(t[0]);
        }
    }
    //D.clear();
}

/**
 * Given two 2d points, calculates the x-coordinate of the intersection point of their lines. Returns
 * -1 if the lines are parallel (or coincident).
 * y = a[0] * x + a[1] * (1 - x);
 * y = b[0] * x + b[1] * (1 - x);
 */
inline
Coord inters_xcoord(Tuple const& a, Tuple const& b){
    Coord shift = a[1] - a[0] - b[1] + b[0];
    if(shift == 0) return -1;
    else{
        Coord diff = a[1] - b[1];
        return diff / shift;
    }
}

struct compare_inters_by_xcoord{
    bool operator() (InterS const& a, InterS const& b) const{
        if(a.xc == b.xc){
            if(min(a.i, a.j) == min(b.i, b.j)){
                return max(a.i, a.j) < max(b.i, b.j);
            }
            return min(a.i, a.j) < min(b.i, b.j);
        }
        return a.xc < b.xc;
    }
};

void get_intersq(DataSet const& SD, DataSet const& kSD, 
    InterSQ & ISQ){

    for(int i = 0; i < SD.size(); i++){
        for(int j = i + 1; j < SD.size(); j++){
            Coord xc = inters_xcoord(SD[i], SD[j]);
            if(xc > 1 || xc < 0){
                cout << "\nerror 1!\n";
                abort();
            }
            ISQ.push_back({xc,i,j});
        }
    }
    for(int i = 0; i < SD.size(); i++){
        for(Tuple const& t: kSD){ //t is a non-skyline
            Coord xc = inters_xcoord(SD[i], t);
            if(0 < xc && xc < 1){ //xc = 0 denotes y-coordinate equal, xc = 1 denotes x-coordinate qual 
                if(SD[i][0] < t[0]) ISQ.push_back({xc,i,-1});
                else if(SD[i][0] > t[0] && SD[i][1] < t[1]) ISQ.push_back({xc,-1,i});
                //else{
                    //cout << "\nerror 2! "<< xc <<"\n";
                    //abort();
                //} 
            }
        }
    }
    sort(ISQ.begin(), ISQ.end(), compare_inters_by_xcoord()); //sort intersections by the x-coordinate in the acscending order
    //kSD.clear();
}

void init_matrix(int const& s, int const& r, vector<vector<Sol>> &MAT){
    for(int i = 0; i < s; i++){
        vector<Sol> row;
        for(int j = 0; j < r; j++){
            vector<SID> sids(1, i); //initial length 1 vector
            row.push_back({sids,0});
        }
        MAT.push_back(row);
    }
}

Sol process_intersq(int const& s, int const& r, int const& k, 
    vector<vector<Sol>> & MAT, SkyRankMap& SRM, InterSQ const&ISQ){
    vector<Coord> lastx(s, 0); 
    //from i to the x-coordinate of the last intersection about i-th skyline tuple
    for(InterS is: ISQ){
        if(is.j == -1){
            if(SRM[is.i] <= k){
                HP delta = get_penalty(SRM[is.i],k) * get_cumul_pr(lastx[is.i], is.xc);
                for(int h = 0; h < r; h++) MAT[is.i][h].hp += delta;
            }
            lastx[is.i] = is.xc;
            SRM[is.i] += 1;
        }else if(is.i == -1){
            if(SRM[is.j] <= k){
                HP delta = get_penalty(SRM[is.j],k) * get_cumul_pr(lastx[is.j], is.xc);
                for(int h = 0; h < r; h++) MAT[is.j][h].hp += delta;
            }
            lastx[is.j] = is.xc;
            SRM[is.j] -= 1;
        }else{
            if(SRM[is.i] + 1 != SRM[is.j]){
                cout << "\nerror 3!\n";
                abort();
            }
            if(SRM[is.i] <= k){
                HP delta = get_penalty(SRM[is.i],k) * get_cumul_pr(lastx[is.i], is.xc);
                for(int h = 0; h < r; h++) MAT[is.i][h].hp += delta;
            }
            if(SRM[is.j] <= k){
                HP delta = get_penalty(SRM[is.j],k) * get_cumul_pr(lastx[is.j], is.xc);
                for(int h = 0; h < r; h++) MAT[is.j][h].hp += delta;
            }

            for(int h = 1; h < r; h++){
                if(MAT[is.j][h].hp < MAT[is.i][h-1].hp){
                    MAT[is.j][h].sids.assign(MAT[is.i][h-1].sids.begin(), MAT[is.i][h-1].sids.end());
                    MAT[is.j][h].sids.push_back(is.j);
                    MAT[is.j][h].hp = MAT[is.i][h-1].hp;
                }
            }

            lastx[is.i] = is.xc, lastx[is.j] = is.xc;
            SRM[is.i] += 1, SRM[is.j] -= 1;
        }
    }
    for(int i = 0; i < s; i++){
        HP delta = get_penalty(SRM[i],k) * get_cumul_pr(lastx[i], 1);
        MAT[i][r-1].hp += delta;
    }
    HP mhp = 0;
    SID rid = -1;
    for(int i = 0; i < s; i++){
        if(MAT[i][r-1].hp > mhp){
            mhp = MAT[i][r-1].hp;
            rid = i;
        }
    }
    return MAT[rid][r-1];
}