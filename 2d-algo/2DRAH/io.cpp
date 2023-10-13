#include "struct.hpp"

void read_file(DataSet & D, int const& n, string const& fileName){
	FILE* fp;
	if((fp = fopen(fileName.c_str(), "r")) == NULL){
		cout<< "Cannot open the data file %s.\n";
		exit(0);
	}

    int nmax,dmax;
    fscanf(fp, "%d", &nmax);
    fscanf(fp, "%d", &dmax);
    //cout << "nmax: " << nmax << ", dmax: " << dmax <<"\n";

    if(nmax < n || dmax < 2){
		cout<< "Input Parameter Error\n";
		exit(0);
	}

	for(int i = 0; i < n; i++){
        Tuple t;
        for(int j = 0; j < dmax; j++){
            double temp;
            fscanf(fp, "%lf", &temp);
            if(j < 2) t[j] = temp;
        }
        D.push_back(t);
	}
	fclose(fp);
}
