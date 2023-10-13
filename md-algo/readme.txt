Readme
=========================
This package contains all source codes for 
MDRAH, MDRAH+, MDRKH, HDRRM and MDRMS

Usage Step
==========
a. Compilation
MDRAH:	g++ -std=c++11 -O2 -g mdrah.cpp
MDRAH+:	g++ -std=c++11 -O2 -g mdrahp.cpp
MDRKH:	g++ -std=c++11 -O2 -g mdrkh.cpp
HDRRM:	g++ -std=c++11 -O2 -g hdrrm.cpp
MDRMS:	g++ -std=c++11 -O2 -g mdrms.cpp
b. Execution
1)For MDRAH
	./a.out data n d r alpha dis_flag epsilon
 * data : data name, e.g., "anti3d", "anti4d", "anti5d", "anti6d", "Colors", "NBA", "Household";
 * n : number of tuple
 * d : dimension
 * r : output size
 * alpha : threshold
 * dis_flag : distribution flag, "1" represents a uniform random distribution on \Delta^(d-1), "3" represents an independent normal distribution with mean 0.5 and standard deviation 0.25
 * epsilon : absolute error of hit pro
2)For MDRAH+
	./a.out data n d r alpha dis_flag epsilon gamma
 * data : data name, e.g., "anti3d", "anti4d", "anti5d", "anti6d", "Colors", "NBA", "Household";
 * n : number of tuple
 * d : dimension
 * r : output size
 * alpha : threshold
 * dis_flag : distribution flag, "1" represents a uniform random distribution on \Delta^(d-1), "3" represents an independent normal distribution with mean 0.5 and standard deviation 0.25
 * epsilon : absolute error of hit pro
 * gamma : discretization parameter
3)For MDRKH, HDRRM and MDRMS
	./a.out data n d r alpha dis_flag
 * data : data name, e.g., "anti3d", "anti4d", "anti5d", "anti6d", "Colors", "NBA", "Household";
 * n : number of tuple
 * d : dimension
 * r : output size
 * alpha : threshold
 * dis_flag : distribution flag, "1" represents a uniform random distribution on \Delta^(d-1), "3" represents an independent normal distribution with mean 0.5 and standard deviation 0.25
c. Output
	The output will be shown on the console. 