#pragma once
#include <vector>
#include <iostream>
#include <algorithm> // is_permutation
#include <set>
#include <map>
#include <list>
#include <string>
#include <random>
#include <sys/time.h>
#include <iomanip>
#include <fstream>

#define CONVERTER 1000000

using namespace std;

using Rank = int; //rank regret
using SID = int; //skyline tuple id
using WID = int; //weight vector id
using Coord = double;
using Weight = double;
using Score = double;
using RRatio = double;
using HP = double;

using Tuple = vector<Coord>; 
using DataSet = vector<Tuple>;

using WVector = vector<Weight>;
using Sol = set<SID>;