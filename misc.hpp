//
//  misc.hpp
//  navstokes
//
//  Created by Alexander Seeliger on 10.11.22.
//

#ifndef misc_hpp
#define misc_hpp

#include <stdio.h>
#include <string>

#endif /* misc_hpp */

std::vector<double> read_parameters(std::string fileName);

void write_parameters(std::string fileName,std::vector<double> U, std::vector<double> V, std::vector<double> P, std::vector<unsigned char> FLAG,int imax, int jmax, double xlength, double ylength);

void write_data(std::string fileName,std::vector<double> U, std::vector<double> V, std::vector<double> P, int timesteps);
