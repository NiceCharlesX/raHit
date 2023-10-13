#include "struct.hpp"

void read_file(DataSet & D, int const& n, int const& d, string const& fileName){
    FILE* fp;
	if((fp = fopen(fileName.c_str(), "r")) == NULL){
		cout<< "Cannot open the data file %s.\n";
		exit(0);
	}

    int nmax,dmax;
    fscanf(fp, "%d", &nmax);
    fscanf(fp, "%d", &dmax);

    if(nmax < n || dmax < d){
		cout<< "Input Parameter Error\n";
		exit(0);
	}

	for(int i = 0; i < n; i++){
        Tuple t(d);
        double temp;
        for(int j = 0; j < dmax; j++){
            fscanf(fp, "%lf", &temp);
            if(j < d) t[j] = temp;
        }
        D.push_back(t);
	}
	fclose(fp);
}