#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "misc.hpp"

using namespace std;

int imax = 5;
int jmax = 5;
double xlength = 1.;
double ylength = 1.;
double delx = xlength / imax;
double dely = ylength / jmax;
vector<double> U;
vector<double> V;
vector<double> P;
double t_end = 5.;
double delt = 1.;
double Re;
double tau;
vector<double> F;
vector<double> G;
vector<double> RHS;

int N;
int itermax;
double eps;
double omega,cgamma,beta,Pr,GX,GY,UI,VI,PI;

vector<double> vector_abs(vector<double> vec)
{
    vector<double> abs = vec; // note to myself: changing abs doesnt change vec
    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i] < 0)
            abs[i] *= -1;
    }
    return abs;
}

// void COMP_DELT(double &delt, int imax, int jmax, double delx, double dely, vector<double> U, vector<double> V, double Re, double tau)
// {
//     vector<double> U_abs = vector_abs(U);
//     vector<double> V_abs = vector_abs(V);
//     delt = tau * min({Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2), delx / *max_element(U_abs.begin(), U_abs.end()), dely / *max_element(V_abs.begin(), V_abs.end())});
// }

void SETBCOND(vector<double> &U, vector<double> &V, int imax, int jmax, int wW, int wE, int wN, int wS)
{
    for (int j = 1; j < jmax; j++)
    {
        V[j] = -V[jmax + 2 + j];
        V[(imax + 1) * (jmax + 2) + j] = -V[imax * (jmax + 2) + j];
    }
    for (int i = 1; i < imax; i++)
    {
        U[i * (jmax + 2)] = -U[i * (jmax + 2) + 1];
        // U[i * (jmax + 2) + jmax + 1] = -U[i * (jmax + 2) + jmax];
        U[i * (jmax + 2) + jmax + 1] = 0.5 -U[i * (jmax + 2) + jmax];
    }
}

void COMP_FG(const vector<double> &U, const vector<double> &V, vector<double> &F, vector<double> &G, int imax, int jmax, double delt, double delx, double dely, double gx, double gy, double gamma, double Re)
{
    double d2UdX2, d2UdY2, dU2dX, dUVdY, d2VdX2, d2VdY2, dUVdX, dV2dY;
    for (int i = 1; i <= imax - 1; i++)
        for (int j = 1; j <= jmax; j++)
        {
            dU2dX = ((U[i * (jmax + 2) + j] + U[(i + 1) * (jmax + 2) + j]) * (U[i * (jmax + 2) + j] + U[(i + 1) * (jmax + 2) + j]) -
                     (U[(i - 1) * (jmax + 2) + j] + U[i * (jmax + 2) + j]) * (U[(i - 1) * (jmax + 2) + j] + U[i * (jmax + 2) + j]) +
                     gamma * (fabs(U[i * (jmax + 2) + j] + U[(i + 1) * (jmax + 2) + j]) * (U[i * (jmax + 2) + j] - U[(i + 1) * (jmax + 2) + j]) -
                              fabs(U[(i - 1) * (jmax + 2) + j] + U[i * (jmax + 2) + j]) * (U[(i - 1) * (jmax + 2) + j] - U[i * (jmax + 2) + j]))) /
                    (4.0 * delx);

            dUVdY = ((V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) * (U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j + 1]) -
                     (V[i * (jmax + 2) + j - 1] + V[(i + 1) * (jmax + 2) + j - 1]) * (U[i * (jmax + 2) + j - 1] + U[i * (jmax + 2) + j]) +
                     gamma * (fabs(V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) * (U[i * (jmax + 2) + j] - U[i * (jmax + 2) + j + 1]) -
                              fabs(V[i * (jmax + 2) + j - 1] + V[(i + 1) * (jmax + 2) + j - 1]) * (U[i * (jmax + 2) + j - 1] - U[i * (jmax + 2) + j]))) /
                    (4.0 * dely);

            d2UdX2 = (U[(i + 1) * (jmax + 2) + j] - 2.0 * U[i * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j]) / delx / delx;

            d2UdY2 = (U[i * (jmax + 2) + j + 1] - 2.0 * U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j - 1]) / dely / dely;

            F[i * (jmax + 2) + j] = U[i * (jmax + 2) + j] + delt * ((d2UdX2 + d2UdY2) / Re - dU2dX - dUVdY + gx);
        }

    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax - 1; j++)
        {
            dUVdX = ((U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j + 1]) * (V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) -
                     (U[(i - 1) * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j + 1]) * (V[(i - 1) * (jmax + 2) + j] + V[i * (jmax + 2) + j]) +
                     gamma * (fabs(U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j + 1]) * (V[i * (jmax + 2) + j] - V[(i + 1) * (jmax + 2) + j]) -
                              fabs(U[(i - 1) * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j + 1]) * (V[(i - 1) * (jmax + 2) + j] - V[i * (jmax + 2) + j]))) /
                    (4.0 * delx);

            dV2dY = ((V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j + 1]) * (V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j + 1]) -
                     (V[i * (jmax + 2) + j - 1] + V[i * (jmax + 2) + j]) * (V[i * (jmax + 2) + j - 1] + V[i * (jmax + 2) + j]) -
                     gamma * (fabs(V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j + 1]) * (V[i * (jmax + 2) + j] - V[i * (jmax + 2) + j + 1]) +
                              fabs(V[i * (jmax + 2) + j - 1] + V[i * (jmax + 2) + j]) * (V[i * (jmax + 2) + j - 1] - V[i * (jmax + 2) + j]))) /
                    (4.0 * dely);

            d2VdX2 = (V[(i + 1) * (jmax + 2) + j] - 2.0 * V[i * (jmax + 2) + j] + V[(i - 1) * (jmax + 2) + j]) / delx / delx;

            d2VdY2 = (V[i * (jmax + 2) + j + 1] - 2.0 * V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j - 1]) / dely / dely;

            G[i * (jmax + 2) + j] = V[i * (jmax + 2) + j] + delt * ((d2VdX2 + d2VdY2) / Re - dUVdX - dV2dY + gy);
        }

    for (int j = 1; j <= jmax; j++)
    {
        F[0 * (jmax + 2) + j] = U[0 * (jmax + 2) + j];
        F[imax * (jmax + 2) + j] = U[imax * (jmax + 2) + j];
    }

    for (int i = 1; i <= imax; i++)
    {
        G[i * (jmax + 2) + 0] = V[i * (jmax + 2) + 0];
        G[i * (jmax + 2) + jmax] = V[i * (jmax + 2) + jmax];
    }
}

void COMP_RHS(const vector<double> &F, const vector<double> &G, vector<double> &RHS, int imax, int jmax, double delt, double delx, double dely)
{
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
            RHS[i * (jmax + 2) + j] = ((F[i * (jmax + 2) + j] - F[(i - 1) * (jmax + 2) + j]) / delx +
                                       (G[i * (jmax + 2) + j] - G[i * (jmax + 2) + j - 1]) / dely) /
                                      delt;
}

int POISSON(vector<double> &P, vector<double> RHS, int imax, int jmax, double delx, double dely, double eps, int itermax, double omg, double &res)
{
    for (int iter = 1; iter <= itermax; iter++)
    {
        double _dx2 = 1 / (delx * delx);
        double _dy2 = 1 / (dely * dely);

        for (int j = 1; j <= jmax; j++)
        {
            P[j] = P[jmax + 2 + j];
            P[(imax + 1) * (jmax + 2) + j] = P[imax * (jmax + 2) + j];
        }

        for (int i = 1; i <= imax; i++)
        {
            P[i * (jmax + 2)] = P[i * (jmax + 2) + 1];
            P[i * (jmax + 2) + jmax + 1] = P[i * (jmax + 2) + jmax];
        }

        for (int i = 1; i <= imax; i++)
            for (int j = 1; j <= jmax; j++)
            {
                P[i * (jmax + 2) + j] = (1 - omg) * P[i * (jmax + 2) + j] +
                                        omg *
                                            ((P[(i + 1) * (jmax + 2) + j] + P[(i - 1) * (jmax + 2) + j]) * _dx2 +
                                             (P[i * (jmax + 2) + j + 1] + P[i * (jmax + 2) + j - 1]) * _dy2 -
                                             RHS[i * (jmax + 2) + j]) /
                                            (2 * (_dx2 + _dy2));
            }

        res = 0;
        double sum = 0;
        for (int i = 1; i <= imax; i++)
            for (int j = 1; j <= jmax; j++)
            {
                sum = (P[(i + 1) * (jmax + 2) + j] - P[i * (jmax + 2) + j] -
                       P[i * (jmax + 2) + j] + P[(i - 1) * (jmax + 2) + j]) *
                          _dx2 +
                      (P[i * (jmax + 2) + j + 1] - P[i * (jmax + 2) + j] -
                       P[i * (jmax + 2) + j] - P[i * (jmax + 2) + j - 1]) *
                          _dy2 -
                      RHS[i * (jmax + 2) + j];
                res += sum * sum;
            }
        res = sqrt(res / (imax * jmax));
        if (res < eps)
            return iter;
    }
    return itermax;
}

void ADAP_UV(vector<double> &U, vector<double> V, const vector<double> &F, const vector<double> &G, const vector<double> &P, int imax, int jmax, double delt, double delx, double dely)
{
    for (int i = 1; i <= imax - 1; i++)
        for (int j = 1; j <= jmax; j++)
        {
            U[i * (jmax + 2) + j] = F[i * (jmax + 2) + j] - (P[(i + 1) * (jmax + 2) + j] - P[i * (jmax + 2) + j]) * delt / delx;
        }
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax - 1; j++)
        {
            V[i * (jmax + 2) + j] = G[i * (jmax + 2) + j] - (P[i * (jmax + 2) + j + 1] - P[i * (jmax + 2) + j]) * delt / dely;
        }
}

void set_parameters(string fileName, double &xlength, double &ylength, int &imax, int &jmax, double &delx, double &dely,
               double &t_end, double &delt, double &tau, int &N,
               int &itermax, double &eps, double &omega, double &cgamma,
               double &Re,double &Pr,double &beta,double &GX,double &GY,
               double &UI,double &VI,double &PI
){
    double* params = read_parameters(fileName);
    xlength=params[0];
    ylength=params[1];
    imax=int(params[2]);
    jmax=int(params[3]);
    delx=params[4];
    dely=params[5];
    t_end=params[6];
    delt=params[7];
    tau=params[8];
    N=int(params[9]);
    itermax=int(params[10]);
    eps=params[11];
    omega=params[12];
    cgamma=params[13];
    Re=params[14];
    Pr=params[15];
    beta=params[16];
    GX=params[17];
    GY=params[18];
    UI=params[19];
    VI=params[20];
    PI=params[21];


}

int main()
{
    
    set_parameters("parameterfiles/test1.txt",xlength,ylength,imax,jmax,delx,dely,t_end,delt,tau,N,itermax,eps,omega,cgamma,Re,Pr,beta,GX,GY,UI,VI,PI);
    cout << xlength << ylength << imax << jmax << delx << dely << t_end << delt<<tau<<N<<itermax<<eps<<omega<<cgamma<<Re<<Pr<<beta<<GX<<GY<< endl;
    cout << UI<< endl;
    cout << VI<<endl;
    cout << PI<<endl;

    imax=5;
    jmax=6;
    vector<double> U((imax+2)*(jmax+2),UI);
    vector<double> V((imax+2)*(jmax+2),VI);
    vector<double> P((imax+2)*(jmax+2),PI);
    cout << U.size() << V.size() << P.size()<< endl;

    write_parameters("newdatafile.txt",U,V,P,imax,jmax,xlength,ylength);

//     int n = 0;
//    double tend = 10;
//    for (double t = 0; t < tend; t += delt)
//    {
//        COMP_DELT(delt, imax, jmax, delx, dely, U, V, Re, tau);
//        SETBCOND(U, V, imax, jmax,0, 0, 0, 0);
//        COMP_FG(U, V, F, G, imax, jmax, delt, delx, dely, gx, gy, gamma, Re);
//        COMP_RHS(F, G, RHS, imax, jmax, delt, delx, dely);
//        POISSON(P, RHS, imax, jmax, delx, dely, eps, itermax, omg, res);
//        ADAP_UV(U, V, F, G, P, imax, jmax, delt, delx, dely);
//        n += 1;
//    }
//
//    cout << "finished in " << n << " steps\n";
//    for (int x : U)
//        cout << x << " ";
//    cout << "\n";
//    for (int x : V)
//        cout << x << " ";
//    cout << "\n";
//    for (int x : P)
//        cout << x << " ";
//    return 0;
}