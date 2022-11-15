#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
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
vector<double> F;
vector<double> G;

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

void COMP_DELT(double &delt, int imax, int jmax, double delx, double dely, vector<double> U, vector<double> V, double Re, double tau)
{
    vector<double> U_abs = vector_abs(U);
    vector<double> V_abs = vector_abs(V);
    delt = tau * min({Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2), delx / *max_element(U_abs.begin(), U_abs.end()), dely / *max_element(V_abs.begin(), V_abs.end())});
}

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
        U[i * (jmax + 2) + jmax + 1] = -U[i * (jmax + 2) + jmax];
    }
}

void COMP_FG(const vector<double> &U, const vector<double> &V, vector<double> &F, vector<double> &G, int imax, int jmax, double delt, double delx, double dely, double GX, double GY, double gamma, double Re)
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

            dUVdY = ((V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) * (U[i * (jmax + 2) + j] + U[i * (jmax + 2) + (j + 1)]) -
                     (V[i * (jmax + 2) + (j - 1)] + V[(i + 1) * (jmax + 2) + (j - 1)]) * (U[i * (jmax + 2) + (j - 1)] + U[i * (jmax + 2) + j]) +
                     gamma * (fabs(V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) * (U[i * (jmax + 2) + j] - U[i * (jmax + 2) + (j + 1)]) -
                              fabs(V[i * (jmax + 2) + (j - 1)] + V[(i + 1) * (jmax + 2) + (j - 1)]) * (U[i * (jmax + 2) + (j - 1)] - U[i * (jmax + 2) + j]))) /
                    (4.0 * dely);

            d2UdX2 = (U[(i + 1) * (jmax + 2) + j] - 2.0 * U[i * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j]) / delx / delx;

            d2UdY2 = (U[i * (jmax + 2) + (j + 1)] - 2.0 * U[i * (jmax + 2) + j] + U[i * (jmax + 2) + (j - 1)]) / dely / dely;

            F[i * (jmax + 2) + j] = U[i * (jmax + 2) + j] + delt * ((d2UdX2 + d2UdY2) / Re - dU2dX - dUVdY + GX)
        }

    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax - 1; j++)
        {
            dUVdX = ((U[i * (jmax + 2) + j] + U[i * (jmax + 2) + (j + 1)]) * (V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) -
                     (U[(i - 1) * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j + 1]) * (V[(i - 1) * (jmax + 2) + j] + V[i * (jmax + 2) + j]) +
                     gamma * (fabs(U[i * (jmax + 2) + j] + U[i * (jmax + 2) + (j + 1)]) * (V[i * (jmax + 2) + j] - V[(i + 1) * (jmax + 2) + j]) -
                              fabs(U[(i - 1) * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j + 1]) * (V[(i - 1) * (jmax + 2) + j] - V[i * (jmax + 2) + j]))) /
                    (4.0 * delx);

            dV2dY = ((V[i * (jmax + 2) + j] + V[i * (jmax + 2) + (j + 1)]) * (V[i * (jmax + 2) + j] + V[i * (jmax + 2) + (j + 1)]) -
                     (V[i * (jmax + 2) + (j - 1)] + V[i * (jmax + 2) + j]) * (V[i * (jmax + 2) + (j - 1)] + V[i * (jmax + 2) + j]) -
                     gamma * (fabs(V[i * (jmax + 2) + j] + V[i * (jmax + 2) + (j + 1)]) * (V[i * (jmax + 2) + j] - V[i * (jmax + 2) + (j + 1)]) +
                              fabs(V[i * (jmax + 2) + (j - 1)] + V[i * (jmax + 2) + j]) * (V[i * (jmax + 2) + (j - 1)] - V[i * (jmax + 2) + j]))) /
                    (4.0 * dely);

            d2VdX2 = (V[(i + 1) * (jmax + 2) + j] - 2.0 * V[i * (jmax + 2) + j] + V[(i - 1) * (jmax + 2) + j]) / delx / delx;

            d2VdY2 = (V[i * (jmax + 2) + (j + 1)] - 2.0 * V[i * (jmax + 2) + j] + V[i * (jmax + 2) + (j - 1)]) / dely / dely;

            G[i * (jmax + 2) + j] = V[i * (jmax + 2) + j] + delt * ((d2VdX2 + d2VdY2) / Re - dUVdX - dV2dY + GY);
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

        // help 2

        return 0;
    }