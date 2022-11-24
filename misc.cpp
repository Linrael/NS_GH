#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cctype>
#include <algorithm>
#include "misc.hpp"

using namespace std;

double* read_parameters(string fileName){

    string paramstrings[] = {"xlength","ylength","imax","jmax","delx","dely","t_end","delt","tau","N","itermax","eps","omega","gamma",
                   "Re","Pr","beta","GX","GY","UI","VI","PI"};
    int siz = sizeof(paramstrings)/sizeof(paramstrings[0]);
    cout << siz << endl;
    double params[siz+1]; //I do not know why this +1 (+2,3 etc. is also possible, but +0 is not) is necessary, but without it, 
                        // when accessing the array values in the main file in the function "set_parameters",
                        // the last TWO variables VI and PI are not set correctly but instead are some random number close to 0
                        // like 2.21297e-314.

   ifstream paramfile;
   paramfile.open(fileName);
   if(paramfile.is_open()){

       cout << "test opened" << endl;

       string searchstring;
       string curr_line;
       size_t length_searchstring;
       size_t pos;

       for(int i=0;i<siz;i++){

            paramfile.seekg(0, ios::beg);
            searchstring=paramstrings[i];
            length_searchstring=int(searchstring.length());

            while ( getline(paramfile, curr_line).good () ){
                if ( (pos = curr_line.find (searchstring, 0)) != string::npos ){
                    curr_line=curr_line.erase(0,length_searchstring);
                    curr_line.erase(remove_if(curr_line.begin(), curr_line.end(), ::isspace), curr_line.end());
                    cout << searchstring << ": " << curr_line << endl;
                    params[i]=stod(curr_line);
                    cout << params[i] << endl;
                    break;
                }
            }
        }
    }

   return params;
};

void write_parameters(std::string fileName, vector<double> U, vector<double> V, vector<double> P, int imax, int jmax){

    ofstream outfile;
    outfile.open(fileName);

    if(!outfile){
        cout << "no file to open";
        return;
    } 

    outfile << imax << endl;
    outfile << jmax << endl;

    for(int i=0;i<U.size();i++){
       outfile << U[i] << "/";
    }
    outfile<< endl;

    for(int i=0;i<V.size();i++){
       outfile << V[i] << "/";
    }
    outfile<< endl;

    for(int i=0;i<P.size();i++){
       outfile << P[i] << "/";
    }

}

