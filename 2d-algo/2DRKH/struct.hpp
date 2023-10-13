#pragma once
#include <array>
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <sys/time.h>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <math.h> 

#define CONVERTER 1000000
using namespace std;

using Rank = int; //positive integer, in {1,...,n}
using SID = int; //-1 denotes non-skyline, 0 is the first skyline tuple
using SkyRankMap = vector<Rank>; //from skyline id to rank now

using Coord = double; //in [0,1]
using Tuple = array<Coord, 2>; //two-dimensional tuple
using DataSet = vector<Tuple>;
using HP = double; //hitting probability, in [0,1]
using Score = double;

struct Sol{
    vector<SID> sids; //the skyline tuples in the solution
    HP hp;
};

struct InterS{
    Coord xc; //x-coordinate
    SID i; //-1 denotes non-skyline
    SID j; //-1 denotes non-skyline
};
using InterSQ = vector<InterS>;
