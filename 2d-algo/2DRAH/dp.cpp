#include "struct.hpp"

void search_restricted_DP(int i, int j, Sol & sol, 
    vector<int> const& PhiI, IntervalList const& IL, vector<vector<bool>> const& DDP, vector<vector<bool>> const& DRDP);

void search_DP(int i, int j, Sol & sol, 
    vector<int> const& PhiI, IntervalList const& IL, vector<vector<bool>> const& DDP, vector<vector<bool>> const& DRDP){
    if(j == 0) return;
    //In the following, j > 0
    if(i < j){
        for(int ii = 0; ii <= i; ii++) sol.sids.push_back(IL[ii].sid); 
        return;
    }
    //In the following, i >= j > 0
    if(DDP[i][j] == true) search_restricted_DP(i, j, sol, PhiI, IL, DDP, DRDP);
    else search_DP(i-1, j, sol, PhiI, IL, DDP, DRDP);
}

void search_restricted_DP(int i, int j, Sol & sol, 
    vector<int> const& PhiI, IntervalList const& IL, vector<vector<bool>> const& DDP, vector<vector<bool>> const& DRDP){
    if(j == 1){
        sol.sids.push_back(IL[i].sid);
        return;
    }
    //In the following, j > 1
    if(i < j){
        for(int ii = 0; ii <= i; ii++) sol.sids.push_back(IL[ii].sid);
        return;
    }
    //In the following, i >= j > 1
    sol.sids.push_back(IL[i].sid);
    if(DRDP[i][j] == true) search_restricted_DP(PhiI[i], j-1, sol, PhiI, IL, DDP, DRDP);
    else search_DP(PhiI[i]-1, j-1, sol, PhiI, IL, DDP, DRDP);
}

//Assumption: IL.size() > r >= 1
Sol execute_dynamic_programming(IntervalList const& IL, int const& r){
    vector<int> PhiI(IL.size(), -1); 
    
    int phii = IL.size() - 2; 
    
    for(int i = IL.size() - 1; i > 0; i--){ 
        while(phii > 0 && IL[phii-1].q >= IL[i].p) phii --; 
        PhiI[i] = phii;
    }
    
    vector<vector<bool>> DDP(IL.size(), vector<bool>(r + 1, false)); //direction of DP
    //DDP[i][j] is true: contain i-th interval
    vector<vector<HP>> DP(IL.size(), vector<HP>(r + 1, 0));
    //DDP[i][j] is true: DP[i][j] = RDP[i][j]
    //DDP[i][j] is false: DP[i][j] = DP[i-1][j]

    //restricted DP: must contain i-th interval
    vector<vector<bool>> DRDP(IL.size(), vector<bool>(r + 1, false)); //direction of restricted DP
    //DRDP[i][j] is true: contain an interval intersecting i-th interval, w.l.o.g., contain PhiI[i]-th interval 
    vector<vector<HP>> RDP(IL.size(), vector<HP>(r + 1, 0));

    for(int i = 0; i < IL.size(); i++) RDP[i][1] = IL[i].q - IL[i].p;

    for(int j = 1; j <= r; j++)
    for(int i = 0; i < j; i++)
    {
        DDP[i][j] = DRDP[i][j] = true;
        DP[i][j] = RDP[i][j] = IL[i].q;
    }
    
    for(int j = 1; j <= r; j++)
    for(int i = j; i < IL.size(); i++)
    {
        if(j > 1){
            if(PhiI[i]-1 < 0 || IL[i].q - IL[PhiI[i]].q + RDP[PhiI[i]][j-1] 
                >= IL[i].q - IL[i].p + DP[PhiI[i]-1][j-1]){
                DRDP[i][j] = true;
                RDP[i][j] = IL[i].q - IL[PhiI[i]].q + RDP[PhiI[i]][j-1];
            }else{
                DRDP[i][j] = false;
                RDP[i][j] = IL[i].q - IL[i].p + DP[PhiI[i]-1][j-1];
            }
        }
        if(RDP[i][j] >= DP[i-1][j]){
            DDP[i][j] = true;
            DP[i][j] = RDP[i][j];
        }else{
            DDP[i][j] = false;
            DP[i][j] = DP[i-1][j];
        }
    }

    Sol sol;
    sol.hp = DP[IL.size()-1][r]; 
    search_DP(IL.size()-1, r, sol, PhiI, IL, DDP, DRDP);
    return sol;
}