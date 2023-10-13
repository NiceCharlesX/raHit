#pragma once
#include <array>
#include <vector>
#include <set>
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

using SID = int; //-1 denotes non-skyline, 0 is the first skyline tuple
using CID = int; //-1 denotes non-convex hull, 0 is the first convex hull tuple
using SkyConvMap = vector<CID>; //from SID to CID
using ConvSkyMap = vector<SID>; //from CID to SID

using Coord = double; //in [0,1]
using Tuple = array<Coord, 2>; //two-dimensional tuple
using DataSet = vector<Tuple>;
using HP = double; //alpha-hit probability, in [0,1]
using Score = double;

struct Sol{
    vector<SID> sids; //the skyline tuples in the solution
    HP hp;
};

struct Interval{ //high scoring space in 2D
    SID sid; //corresponding sid
    double p; //left border
    double q; //right border
};
using IntervalList = vector<Interval>;
