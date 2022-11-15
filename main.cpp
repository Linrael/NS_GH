#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

int imax = 100;
int jmax = 100;
double a = 1;
double b = 1;
double delx = a / imax;
double dely = b / jmax;
vector<double> U;
vector<double> V;
double delt = 1.;
double Re;
double tau;

vector<double> vector_abs(vector<double> vec){
    vector<double> abs = vec; // note to myself: changing abs doesnt change vec
    for(int i = 0; i < vec.size(); i++)
    {
        if(vec[i] < 0) abs[i] *= -1;
    }
    return abs;
}

void COMP_DELT(double &delt,int imax,int jmax,double delx,double dely,vector<double> U,vector<double> V,double Re,double tau){
    vector<double> U_abs = vector_abs(U);
    vector<double> V_abs = vector_abs(V);
    delt = tau * min({Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2), delx / *max_element(U_abs.begin(), U_abs.end()), dely / *max_element(V_abs.begin(), V_abs.end())});
}

void SETBCOND(vector<double> U,vector<double> V,int imax,int jmax,int wW,int wE,int wN,int wS){
    for(int j=1; j<jmax; j++){
        V[j] = -V[jmax+2 + j];
        V[(imax + 1)*(jmax+2) + j] = -V[imax*(jmax+2) + j];
    }
    for(int i=1; i<imax; i++){
        U[i*(jmax+2)] = -U[i*(jmax+2) + 1];
        U[i*(jmax+2) + jmax + 1] = -U[i*(jmax + 2) + jmax];
    }
}

int c=3;

void test(int &c){
    c=2;
}

int main()
{
    // vector<double> a = {1,-22,2,-3,5};
    // vector<double> b = vector_abs(a);
    // b.push_back(8);

    // for (double i: b)
    //     cout << i << ' ';

    // cout << 5;

    cout << c;
    test(c);
    cout << c;
    
    return 0;
}