#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "misc.hpp"
#include <math.h>
#include <chrono>
#include <ctime>

// some constant values
const unsigned char C_B = 0b0000; /* This cell is an obstacle/boundary cell */
const unsigned char B_N = 0b0001; /* This obstacle cell has a fluid cell to the north */
const unsigned char B_S = 0b0010; /* This obstacle cell has a fluid cell to the south */
const unsigned char B_W = 0b0100; /* This obstacle cell has a fluid cell to the west */
const unsigned char B_E = 0b1000; /* This obstacle cell has a fluid cell to the east */
const unsigned char B_NW = (B_N | B_W);
const unsigned char B_SW = (B_S | B_W);
const unsigned char B_NE = (B_N | B_E);
const unsigned char B_SE = (B_S | B_E);
const unsigned char B_NSEW = (B_N | B_S | B_E | B_W);

const unsigned char C_F = 0b10000; /* This cell is a fluid cell */

using namespace std;

vector<double> vector_abs(vector<double> vec) {
    vector<double> abs = vec;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] < 0)
            abs[i] *= -1;
    }
    return abs;
}

void
COMP_DELT(double &delt, int imax, int jmax, double delx, double dely, vector<double> &U, vector<double> &V, double Re,
          double tau, int problem) {

    if(problem==3){
        delt=0.001;
        return;
    }
    vector<double> U_abs = vector_abs(U);
    vector<double> V_abs = vector_abs(V);
    double u = *max_element(U_abs.begin(), U_abs.end());
    double v = *max_element(V_abs.begin(), V_abs.end());
    if (u <= 1.e-10) {
        if (v <= 1.e-10) {
            delt = tau * Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2);
        } else {
            delt = tau * min(Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2), dely / v);
        }
    } else {
        if (v <= 1.e-10) {
            delt = tau * min(Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2), delx / u);
        } else {
            delt = min(delx / u, dely / v);
            delt = tau * min(Re / ((1 / (delx * delx) + 1 / (dely * dely)) * 2), delt);
        }
    }


}

void SETBCOND(vector<double> &U, vector<double> &V, int imax, int jmax, int wW, int wE, int wN, int wS) {
    for (int j = 1; j < jmax; j++) {
        V[j] = -V[jmax + 2 + j];
        V[(imax + 1) * (jmax + 2) + j] = -V[imax * (jmax + 2) + j];
    }
    for (int i = 1; i < imax; i++) {
        U[i * (jmax + 2)] = -U[i * (jmax + 2) + 1];
        // U[i * (jmax + 2) + jmax + 1] = -U[i * (jmax + 2) + jmax];
        U[i * (jmax + 2) + jmax + 1] = 0.5 - U[i * (jmax + 2) + jmax];
    }
}

void INITFLAG(vector<unsigned char> &flag, int imax, int jmax) {
    /* Mark the north & south boundary cells */
    for (int i = 0; i <= imax + 1; i++) {
        flag[i * (jmax + 2) + 0] = 0;
        flag[i * (jmax + 2) + jmax + 1] = 0;
    }
    /* Mark the east and west boundary cells */
    for (int j = 1; j <= jmax; j++) {
        flag[0 * (jmax + 2) + j] = 0;
        flag[(imax + 1) * (jmax + 2) + j] = 0;
    }

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            if (!(flag[i * (jmax + 2) + j] & C_F)) {
                if (flag[(i - 1) * (jmax + 2) + j] & C_F)
                    flag[i * (jmax + 2) + j] |= B_W;
                if (flag[(i + 1) * (jmax + 2) + j] & C_F)
                    flag[i * (jmax + 2) + j] |= B_E;
                if (flag[i * (jmax + 2) + j - 1] & C_F)
                    flag[i * (jmax + 2) + j] |= B_S;
                if (flag[i * (jmax + 2) + j + 1] & C_F)
                    flag[i * (jmax + 2) + j] |= B_N;
            }
        }
    }
}

void SETBCONDnew(vector<double> &U, vector<double> &V, vector<double> &P, vector<unsigned char> &FLAG, int imax, int jmax,
            int wW, int wE, int wN, int wS, int problem) {
    /* The flags wW,wE,wN, and wS can have the values:
    //  1 = slip               2 = no-slip
    //  3 = outflow            4 = periodic
    //  Moreover, no-slip conditions are set at internal obstacle cells
    //  by default.
    */
    if (problem == 1) {
        //inflow from the left for upper half
        for (int j = int(0.5 * jmax); j <= jmax + 1; j++) {
            U[1 * (jmax + 2) + j] = 1.;
        }
    }

    if (problem == 4) {
        //inflow from the left
        for (int j = floor(jmax/5); j <= jmax - floor(3*jmax/5); j++) {
            U[1 * (jmax + 2) + j] = 1.;
        }
    }


    switch (wW) {
        case 1:
            for (int j = 0; j <= jmax + 1; j++) {
                U[0 * (jmax + 2) + j] = 0.0;
                V[0 * (jmax + 2) + j] = V[1 * (jmax + 2) + j];
            }
            break;

        case 2:
            for (int j = 0; j <= jmax + 1; j++) {
                U[0 * (jmax + 2) + j] = 0.0;
                V[0 * (jmax + 2) + j] = (-1.0) * V[1 * (jmax + 2) + j];
            }
            break;

        case 3:
            for (int j = 0; j <= jmax + 1; j++) {
                U[0 * (jmax + 2) + j] = U[1 * (jmax + 2) + j];
                V[0 * (jmax + 2) + j] = V[1 * (jmax + 2) + j];
            }
            break;

        case 4:
            for (int j = 0; j <= jmax + 1; j++) {
                U[0 * (jmax + 2) + j] = U[(imax - 1) * (jmax + 2) + j];
                V[0 * (jmax + 2) + j] = V[(imax - 1) * (jmax + 2) + j];
                V[1 * (jmax + 2) + j] = V[imax * (jmax + 2) + j];
                P[1 * (jmax + 2) + j] = P[imax * (jmax + 2) + j];
            }
            break;
    }

    switch (wE) {
        case 1:
            for (int j = 0; j <= jmax + 1; j++) {
                U[imax * (jmax + 2) + j] = 0.0;
                V[(imax + 1) * (jmax + 2) + j] = V[imax * (jmax + 2) + j];
            }
            break;

        case 2:
            for (int j = 0; j <= jmax + 1; j++) {
                U[imax * (jmax + 2) + j] = 0.0;
                V[(imax + 1) * (jmax + 2) + j] = (-1.0) * V[imax * (jmax + 2) + j];
            }
            break;

        case 3:
            for (int j = 0; j <= jmax + 1; j++) {
                U[imax * (jmax + 2) + j] = U[(imax - 1) * (jmax + 2) + j];
                V[(imax + 1) * (jmax + 2) + j] = V[imax * (jmax + 2) + j];
            }
            break;

        case 4:
            for (int j = 0; j <= jmax + 1; j++) {
                U[imax * (jmax + 2) + j] = U[1 * (jmax + 2) + j];
                V[(imax + 1) * (jmax + 2) + j] = V[2 * (jmax + 2) + j];
            }
            break;
    }

    switch (wN) {
        case 1:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + jmax] = 0.0;
                U[i * (jmax + 2) + jmax + 1] = U[i * (jmax + 2) + jmax];
            }
            break;

        case 2:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + jmax] = 0.0;
                U[i * (jmax + 2) + jmax + 1] = (-1.0) * U[i * (jmax + 2) + jmax];
            }
            break;

        case 3:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + jmax] = V[i * (jmax + 2) + jmax - 1];
                U[i * (jmax + 2) + jmax + 1] = U[i * (jmax + 2) + jmax];
            }
            break;

        case 4:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + jmax] = V[i * (jmax + 2) + 1];
                U[i * (jmax + 2) + jmax + 1] = U[i * (jmax + 2) + 2];
            }
            break;
    }

    switch (wS) {
        case 1:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + 0] = 0.0;
                U[i * (jmax + 2) + 0] = U[i * (jmax + 2) + 1];
            }
            break;

        case 2:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + 0] = 0.0;
                U[i * (jmax + 2) + 0] = (-1.0) * U[i * (jmax + 2) + 1];
            }
            break;

        case 3:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + 0] = V[i * (jmax + 2) + 1];
                U[i * (jmax + 2) + 0] = U[i * (jmax + 2) + 1];
            }
            break;

        case 4:
            for (int i = 0; i <= imax + 1; i++) {
                V[i * (jmax + 2) + 0] = V[i * (jmax + 2) + jmax - 1];
                U[i * (jmax + 2) + 0] = U[i * (jmax + 2) + jmax - 1];
                U[i * (jmax + 2) + 1] = U[i * (jmax + 2) + jmax];
                P[i * (jmax + 2) + 1] = P[i * (jmax + 2) + jmax];
            }
            break;
    }

    if(problem == 3){
        for (int i = 0; i <= imax+1; ++i){
            U[i*(jmax+2)+jmax+1]=1-U[i*(jmax+2)+jmax];
        }//U[i * (jmax + 2) + jmax + 1] = (-1.0) * U[i * (jmax + 2) + jmax];
    };



    // Now apply BC for inner cells: (only noslip conditions for inner cells)
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
            if (FLAG[i * (jmax + 2) + j] & 0b01111) {
                //
                //  The mask 0x000f filters the obstacle cells adjacent to fluid cells.
                //
                switch (FLAG[i * (jmax + 2) + j]) {
                    case B_N: {
                        V[i * (jmax + 2) + j] = 0.0;
                        U[i * (jmax + 2) + j] = -U[i * (jmax + 2) + j + 1];
                        U[(i - 1) * (jmax + 2) + j] = -U[(i - 1) * (jmax + 2) + j + 1];
                        break;
                    }
                    case B_E: {
                        U[i * (jmax + 2) + j] = 0.0;
                        V[i * (jmax + 2) + j] = -V[(i + 1) * (jmax + 2) + j];
                        V[i * (jmax + 2) + j - 1] = -V[(i + 1) * (jmax + 2) + j - 1];
                        break;
                    }
                    case B_S: {
                        V[i * (jmax + 2) + j - 1] = 0.0;
                        U[i * (jmax + 2) + j] = -U[i * (jmax + 2) + j - 1];
                        U[(i - 1) * (jmax + 2) + j] = -U[(i - 1) * (jmax + 2) + j - 1];
                        break;
                    }
                    case B_W: {
                        U[(i - 1) * (jmax + 2) + j] = 0.0;
                        V[i * (jmax + 2) + j] = -V[(i - 1) * (jmax + 2) + j];
                        V[i * (jmax + 2) + j - 1] = -V[(i - 1) * (jmax + 2) + j - 1];
                        break;
                    }
                    case B_NE: {
                        V[i * (jmax + 2) + j] = 0.0;
                        U[i * (jmax + 2) + j] = 0.0;
                        V[i * (jmax + 2) + j - 1] = -V[(i + 1) * (jmax + 2) + -1];
                        U[(i - 1) * (jmax + 2) + j] = -U[(i - 1) * (jmax + 2) + j + 1];
                        break;
                    }
                    case B_SE: {
                        V[i * (jmax + 2) + j - 1] = 0.0;
                        U[i * (jmax + 2) + j] = 0.0;
                        V[i * (jmax + 2) + j] = -V[(i + 1) * (jmax + 2) + j];
                        U[(i - 1) * (jmax + 2) + j] = -U[(i - 1) * (jmax + 2) + j - 1];
                        break;
                    }
                    case B_SW: {
                        V[i * (jmax + 2) + j - 1] = 0.0;
                        U[(i - 1) * (jmax + 2) + j] = 0.0;
                        V[i * (jmax + 2) + j] = -V[(i - 1) * (jmax + 2) + j];
                        U[i * (jmax + 2) + j] = -U[i * (jmax + 2) + j - 1];
                        break;
                    }
                    case B_NW: {
                        V[i * (jmax + 2) + j] = 0.0;
                        U[(i - 1) * (jmax + 2) + j] = 0.0;
                        V[i * (jmax + 2) + j - 1] = -V[(i - 1) * (jmax + 2) + j - 1];
                        U[i * (jmax + 2) + j] = -U[i * (jmax + 2) + j + 1];
                        break;
                    }
                }
            }
}

void COMP_FG(const vector<double> &U, const vector<double> &V, vector<double> &F, vector<double> &G,
             vector<unsigned char> &FLAG, int imax, int jmax, double delt, double delx, double dely, double gx,
             double gy, double cgamma, double Re) {
    // only if both adjacent cells are fluid
    double d2UdX2, d2UdY2, dU2dX, dUVdY, d2VdX2, d2VdY2, dUVdX, dV2dY;
    for (int i = 1; i <= imax - 1; i++)
        for (int j = 1; j <= jmax; j++) {
            if ((FLAG[i * (jmax + 2) + j] & C_F) && (FLAG[(i + 1) * (jmax + 2) + j] & C_F)) {
                dU2dX = ((U[i * (jmax + 2) + j] + U[(i + 1) * (jmax + 2) + j]) *
                         (U[i * (jmax + 2) + j] + U[(i + 1) * (jmax + 2) + j]) -
                         (U[(i - 1) * (jmax + 2) + j] + U[i * (jmax + 2) + j]) *
                         (U[(i - 1) * (jmax + 2) + j] + U[i * (jmax + 2) + j]) +
                         cgamma * (fabs(U[i * (jmax + 2) + j] + U[(i + 1) * (jmax + 2) + j]) *
                                   (U[i * (jmax + 2) + j] - U[(i + 1) * (jmax + 2) + j]) -
                                   fabs(U[(i - 1) * (jmax + 2) + j] + U[i * (jmax + 2) + j]) *
                                   (U[(i - 1) * (jmax + 2) + j] - U[i * (jmax + 2) + j]))) /
                        (4.0 * delx);

                dUVdY = ((V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) *
                         (U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j + 1]) -
                         (V[i * (jmax + 2) + j - 1] + V[(i + 1) * (jmax + 2) + j - 1]) *
                         (U[i * (jmax + 2) + j - 1] + U[i * (jmax + 2) + j]) +
                         cgamma * (fabs(V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) *
                                   (U[i * (jmax + 2) + j] - U[i * (jmax + 2) + j + 1]) -
                                   fabs(V[i * (jmax + 2) + j - 1] + V[(i + 1) * (jmax + 2) + j - 1]) *
                                   (U[i * (jmax + 2) + j - 1] - U[i * (jmax + 2) + j]))) /
                        (4.0 * dely);

                d2UdX2 = (U[(i + 1) * (jmax + 2) + j] - 2.0 * U[i * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j]) /
                         delx / delx;

                d2UdY2 = (U[i * (jmax + 2) + j + 1] - 2.0 * U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j - 1]) / dely /
                         dely;

                F[i * (jmax + 2) + j] = U[i * (jmax + 2) + j] + delt * ((d2UdX2 + d2UdY2) / Re - dU2dX - dUVdY + gx);
            } else
                F[i * (jmax + 2) + j] = U[i * (jmax + 2) + j];
        }

    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax - 1; j++) {
            if ((FLAG[i * (jmax + 2) + j] & C_F) && (FLAG[i * (jmax + 2) + j + 1] & C_F)) {
                dUVdX = ((U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j + 1]) *
                         (V[i * (jmax + 2) + j] + V[(i + 1) * (jmax + 2) + j]) -
                         (U[(i - 1) * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j + 1]) *
                         (V[(i - 1) * (jmax + 2) + j] + V[i * (jmax + 2) + j]) +
                         cgamma * (fabs(U[i * (jmax + 2) + j] + U[i * (jmax + 2) + j + 1]) *
                                   (V[i * (jmax + 2) + j] - V[(i + 1) * (jmax + 2) + j]) -
                                   fabs(U[(i - 1) * (jmax + 2) + j] + U[(i - 1) * (jmax + 2) + j + 1]) *
                                   (V[(i - 1) * (jmax + 2) + j] - V[i * (jmax + 2) + j]))) /
                        (4.0 * delx);

                dV2dY = ((V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j + 1]) *
                         (V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j + 1]) -
                         (V[i * (jmax + 2) + j - 1] + V[i * (jmax + 2) + j]) *
                         (V[i * (jmax + 2) + j - 1] + V[i * (jmax + 2) + j]) -
                         cgamma * (fabs(V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j + 1]) *
                                   (V[i * (jmax + 2) + j] - V[i * (jmax + 2) + j + 1]) +
                                   fabs(V[i * (jmax + 2) + j - 1] + V[i * (jmax + 2) + j]) *
                                   (V[i * (jmax + 2) + j - 1] - V[i * (jmax + 2) + j]))) /
                        (4.0 * dely);

                d2VdX2 = (V[(i + 1) * (jmax + 2) + j] - 2.0 * V[i * (jmax + 2) + j] + V[(i - 1) * (jmax + 2) + j]) /
                         delx / delx;

                d2VdY2 = (V[i * (jmax + 2) + j + 1] - 2.0 * V[i * (jmax + 2) + j] + V[i * (jmax + 2) + j - 1]) / dely /
                         dely;

                G[i * (jmax + 2) + j] = V[i * (jmax + 2) + j] + delt * ((d2VdX2 + d2VdY2) / Re - dUVdX - dV2dY + gy);
            } else
                G[i * (jmax + 2) + j] = V[i * (jmax + 2) + j];
        }

    for (int j = 1; j <= jmax; j++) {
        F[0 * (jmax + 2) + j] = U[0 * (jmax + 2) + j];
        F[imax * (jmax + 2) + j] = U[imax * (jmax + 2) + j];
    }

    for (int i = 1; i <= imax; i++) {
        G[i * (jmax + 2) + 0] = V[i * (jmax + 2) + 0];
        G[i * (jmax + 2) + jmax] = V[i * (jmax + 2) + jmax];
    }
}

void COMP_RHS(const vector<double> &F, const vector<double> &G, vector<double> &RHS, vector<unsigned char> &FLAG, int imax,
         int jmax, double delt, double delx, double dely) {
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
            if (FLAG[i * (jmax + 2) + j] & C_F)
                RHS[i * (jmax + 2) + j] = ((F[i * (jmax + 2) + j] - F[(i - 1) * (jmax + 2) + j]) / delx +
                                           (G[i * (jmax + 2) + j] - G[i * (jmax + 2) + j - 1]) / dely) /
                                          delt;
}

int POISSON(vector<double> &P, vector<double> RHS, vector<unsigned char> &FLAG, int imax, int jmax, double delx,
            double dely, double eps, int itermax, double omega, double &res, double t) {
    int iter;
    double _dx2,_dy2;
//    vector<double> atemp(imax * jmax, 0);
    for (iter = 1; iter <= itermax; iter++) {
        _dx2 = 1 / (delx * delx);
        _dy2 = 1 / (dely * dely);

        // external boundaries see eq 3.41
        for (int j = 1; j <= jmax; j++) {
            P[j] = P[jmax + 2 + j];
            P[(imax + 1) * (jmax + 2) + j] = P[imax * (jmax + 2) + j];
        }

        for (int i = 1; i <= imax; i++) {
            P[i * (jmax + 2)] = P[i * (jmax + 2) + 1];
            P[i * (jmax + 2) + jmax + 1] = P[i * (jmax + 2) + jmax];
        }

        // internal boundaries see eq. after 3.53
        for (int i = 1; i <= imax; i++)
            for (int j = 1; j <= jmax; j++)
                if (FLAG[i * (jmax + 2) + j] >= B_N && FLAG[i * (jmax + 2) + j] <= B_SE)
                    switch (FLAG[i * (jmax + 2) + j]) {
                        case B_N: {
                            P[i * (jmax + 2) + j] = P[i * (jmax + 2) + j + 1];
                            break;
                        }
                        case B_E: {
                            P[i * (jmax + 2) + j] = P[(i + 1) * (jmax + 2) + j];
                            break;
                        }
                        case B_S: {
                            P[i * (jmax + 2) + j] = P[i * (jmax + 2) + j - 1];
                            break;
                        }
                        case B_W: {
                            P[i * (jmax + 2) + j] = P[(i - 1) * (jmax + 2) + j];
                            break;
                        }
                        case B_NE: {
                            P[i * (jmax + 2) + j] = 0.5 * (P[i * (jmax + 2) + j + 1] + P[(i + 1) * (jmax + 2) + j]);
                            break;
                        }
                        case B_SE: {
                            P[i * (jmax + 2) + j] = 0.5 * (P[i * (jmax + 2) + j - 1] + P[(i + 1) * (jmax + 2) + j]);
                            break;
                        }
                        case B_SW: {
                            P[i * (jmax + 2) + j] = 0.5 * (P[i * (jmax + 2) + j - 1] + P[(i - 1) * (jmax + 2) + j]);
                            break;
                        }
                        case B_NW: {
                            P[i * (jmax + 2) + j] = 0.5 * (P[i * (jmax + 2) + j + 1] + P[(i - 1) * (jmax + 2) + j]);
                            break;
                        }
                    }

        // for all fluid cells
        //#pragma omp parallel for #TODO
        for (int i = 1; i <= imax; i++)
            for (int j = 1; j <= jmax; j++) {
                if (FLAG[i * (jmax + 2) + j] & C_F) {
                    P[i * (jmax + 2) + j] = (1 - omega) * P[i * (jmax + 2) + j] +
                                            omega *
                                            ((P[(i + 1) * (jmax + 2) + j] + P[(i - 1) * (jmax + 2) + j]) * _dx2 +
                                             (P[i * (jmax + 2) + j + 1] + P[i * (jmax + 2) + j - 1]) * _dy2 -
                                             RHS[i * (jmax + 2) + j]) /
                                            (2 * (_dx2 + _dy2));
                }
            }

        // calculate residual
        res = 0;
        double sum = 0;
        int ij = 0;
        for (int i = 1; i <= imax; i++)
            for (int j = 1; j <= jmax; j++) {
                if (FLAG[i * (jmax + 2) + j] & C_F) {
//                    atemp[(i - 1) * jmax + j] = (P[(i + 1) * (jmax + 2) + j] - P[i * (jmax + 2) + j] -
//                                                 P[i * (jmax + 2) + j] + P[(i - 1) * (jmax + 2) + j]) *
//                                                _dx2 +
//                                                (P[i * (jmax + 2) + j + 1] - P[i * (jmax + 2) + j] -
//                                                 P[i * (jmax + 2) + j] - P[i * (jmax + 2) + j - 1]) *
//                                                _dy2;

                    sum = (P[(i + 1) * (jmax + 2) + j] - P[i * (jmax + 2) + j] -
                        P[i * (jmax + 2) + j] + P[(i - 1) * (jmax + 2) + j]) *
                       _dx2 +
                       (P[i * (jmax + 2) + j + 1] - P[i * (jmax + 2) + j] -
                        P[i * (jmax + 2) + j] + P[i * (jmax + 2) + j - 1]) *
                       _dy2 - RHS[i * (jmax + 2) + j];
                    res += sum * sum;
                    ij++;
                }
            }
        res = sqrt(res / ij);
        if (res < eps) {
            cout << iter << endl;
            return iter;
        }
    }
//    if (t > 0.7) {
//        cout << res << " ; " << iter << endl;
//        int counter = 0;
//        for (int myloop = 0; myloop < imax * jmax; myloop++) {
//            if (myloop % jmax == 0) {
//                counter++;
//                cout << endl;
//                cout << "P" << counter << ": ";
//            }
//            cout << atemp[myloop] << " ";
//        }
//        counter = 0;
//        for (int myloop = 0; myloop < imax * jmax; myloop++) {
//            if (myloop % jmax == 0) {
//                counter++;
//                cout << endl;
//                cout << "rhs" << counter << ": ";
//
//            }
//            cout << btemp[myloop] << " ";
//        }
//    }
    cout << res << " ; " << iter << endl;
    return itermax;
}

void ADAP_UV(vector<double> &U, vector<double> &V, const vector<double> &F, const vector<double> &G, const vector<double> &P,
        const vector<unsigned char> &FLAG, int imax, int jmax, double delt, double delx, double dely) {
    // only if both adjacent cells are fluid
    for (int i = 1; i <= imax - 1; i++)
        for (int j = 1; j <= jmax; j++)
            if ((FLAG[i * (jmax + 2) + j] & C_F) && (FLAG[(i + 1) * (jmax + 2) + j] & C_F))
                U[i * (jmax + 2) + j] =
                        F[i * (jmax + 2) + j] - (P[(i + 1) * (jmax + 2) + j] - P[i * (jmax + 2) + j]) * delt / delx;

    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax - 1; j++)
            if ((FLAG[i * (jmax + 2) + j] & C_F) && (FLAG[(i + 1) * (jmax + 2) + j] & C_F))
                V[i * (jmax + 2) + j] =
                        G[i * (jmax + 2) + j] - (P[i * (jmax + 2) + j + 1] - P[i * (jmax + 2) + j]) * delt / dely;
}

void set_parameters(string fileName, double &xlength, double &ylength, int &imax, int &jmax, double &delx, double &dely,
                    double &t_end, double &delt, double &tau, int &N,
                    int &itermax, double &eps, double &omega, double &cgamma,
                    double &Re, double &Pr, double &beta, double &gx, double &gy,
                    double &UI, double &VI, double &PI) {
    vector<double> params = read_parameters(fileName);
    xlength = params[0];
    ylength = params[1];
    imax = int(params[2]);
    jmax = int(params[3]);
    delx = params[4];
    dely = params[5];
    t_end = params[6];
    delt = params[7];
    tau = params[8];
    N = int(params[9]);
    itermax = int(params[10]);
    eps = params[11];
    omega = params[12];
    cgamma = params[13];
    Re = params[14];
    Pr = params[15];
    beta = params[16];
    gx = params[17];
    gy = params[18];
    UI = params[19];
    VI = params[20];
    PI = params[21];
}

void SIMPLEGEOMETRY(vector<unsigned char> &FLAG, int i1, int i2, int j1, int j2, int imax, int jmax) {
    if (i1 < 1 || j1 < 1 || i2 > imax || j2 > jmax) {
        cout << "FEHLER IN GEOMETRY" << endl;
        return;
    }
    for (int i = i1; i <= i2; i++) {
        for (int j = j1; j <= j2; j++) {
            FLAG[i * (jmax + 2) + j] = C_B;
        }
    }
}

void roehren(vector<unsigned char> &FLAG, int imax, int jmax){
    int inout_i = floor(imax/4);
    int small_j = floor(jmax/5);
    int large_j = floor(3*jmax/5);
    int diameter = floor(min(imax, jmax) / 5);
    
    SIMPLEGEOMETRY(FLAG, 1, inout_i, 1, small_j, imax, jmax);
    SIMPLEGEOMETRY(FLAG, 1, inout_i, jmax-large_j, jmax, imax, jmax);
    SIMPLEGEOMETRY(FLAG, inout_i + diameter, imax - inout_i - diameter, diameter, jmax - diameter, imax, jmax);
    SIMPLEGEOMETRY(FLAG, imax - inout_i, imax, 1, large_j, imax, jmax);
    SIMPLEGEOMETRY(FLAG, imax - inout_i, imax, jmax - small_j, jmax, imax, jmax);
}

vector<int> koutofn(int k, int n) {
    //returns k evenly spaced integers from a sequence of n consecutive integers
    vector<int> ret(k, 0);
    for (int i = 0; i < k; i++) {
        ret[i] = floor(i * n / k) + floor(n / (2 * k));
    }
    return ret;
}

void trichter(vector<unsigned char> &FLAG, int linkegrenze, int rechtegrenze, const string &seite, int imax, int jmax) {
    int differenz = abs(linkegrenze - rechtegrenze);
    if (differenz == 0) {
        cout << "bitte simplegeometry nutzen" << endl;
        return;
    }
    int steigung = int(floor(imax / differenz));
    int zuwenig = imax - steigung * differenz;
    int abstand_zur_seite = abs(jmax - linkegrenze);

    vector<int> schritte_bis_einshoch(differenz, steigung);

    if (zuwenig != 0) {
        vector<int> muss_eins_dazu = koutofn(zuwenig, differenz);
        vector<int> schritte_bis_einshoch(differenz, 0);
        int j = 0;
        for (int i = 0; i < differenz; i++) {
            if (i == muss_eins_dazu[j]) {
                schritte_bis_einshoch[i]++;
                j++;
            }
        }
    }

    int x = 1;
    int x_current = x;
    int y_height = 0;
    int index = 0;
    for (index = 0; index < differenz; index++) {
        for (x_current = x; (x_current < x + schritte_bis_einshoch[index]) && (x_current <= imax); x_current++) {
            if (seite == "oben") {
                for (int j = 0; j <= y_height + abstand_zur_seite; j++) {
                    FLAG[x_current * (jmax + 2) + j + 1] = 0;
                }
            } else if (seite == "unten") {
                for (int j = 0; j <= y_height + (jmax - abstand_zur_seite); j++) {
                    FLAG[(x_current + 1) * (jmax + 2) - j - 1] = 0;
                }
            } else {
                cout << "sinnlos";
            };
        }
        x = x_current;
        y_height++;
    }
}

int main() {
    int imax;
    int jmax;
    double xlength;
    double ylength;
    double delx;
    double dely;
    double t_end;
    double delt;
    double Re;
    double tau;
    char problem = 4;

    int N; // Number of particle lines
    int itermax;
    double eps;
    double omega, cgamma, beta, Pr, gx, gy, UI, VI, PI, res;

    string outputfile;

    set_parameters("parameterfiles/test1.txt", xlength, ylength, imax, jmax, delx, dely, t_end, delt, tau, N,
                   itermax, eps, omega, cgamma, Re, Pr, beta, gx, gy, UI, VI, PI);
    cout << xlength << ylength << imax << jmax << delx << dely << t_end << delt << tau << N << itermax << eps << omega
         << cgamma << Re << Pr << beta << gx << gy << endl;
    cout << UI << endl;
    cout << VI << endl;
    cout << PI << endl;

    delx = xlength / imax;
    dely = ylength / jmax;

    vector<double> F((imax + 2) * (jmax + 2), 0.);
    vector<double> G((imax + 2) * (jmax + 2), 0.);
    vector<double> RHS((imax + 2) * (jmax + 2), 0.);

    vector<double> U((imax + 2) * (jmax + 2), UI);
    vector<double> V((imax + 2) * (jmax + 2), VI);
    vector<double> P((imax + 2) * (jmax + 2), PI);
    vector<unsigned char> FLAG((imax + 2) * (jmax + 2), 0b10000);

    int wW, wE, wN, wS;

    switch (problem) {
        case 0: //square in the middle
            SIMPLEGEOMETRY(FLAG, 40, 50, 40, 50, imax, jmax);
            break;
        case 1: //square bottom left corner
            SIMPLEGEOMETRY(FLAG, 1, floor(0.2 * imax), 1, floor(0.5 * jmax), imax, jmax);
            for (int i = 0; i <= imax + 1; ++i) {
                for (int j = floor(0.5 * jmax); j <= jmax + 1; j++) {
                    U[i * (jmax + 2) + j] = 1.;
                }
            }
            wW=3;
            wE=3;
            wN=2;
            wS=2;
            break;
        case 2: // Trichter
            trichter(FLAG, 90, 70, "oben", imax, jmax);
            trichter(FLAG, 10, 30, "unten", imax, jmax);
            break;
        case 3: //Self Convergence Lid Driven Cavity
            wW=2;
            wE=2;
            wN=2;
            wS=2;
            break;
        case 4: // Röhre
            roehren(FLAG, imax, jmax);
            wW=3;
            wE=3;
            wN=2;
            wS=2;
            break;

        default:
            break;
    }

    INITFLAG(FLAG, imax, jmax);

    write_parameters("finaldata/liddrivencavity200.txt", U, V, P, FLAG, imax, jmax, xlength, ylength, delt);
    int t1 = time(0);

    int n = 0;
    for (double t = 0; t < t_end; t += delt) {
        cout << t << endl;
        /*
        COMP_DELT(delt, imax, jmax, delx, dely, U, V, Re, tau);
        SETBCOND(U, V, imax, jmax, 0, 0, 0, 0);
        COMP_FG(U, V, F, G, imax, jmax, delt, delx, dely, gx, gy, cgamma, Re);
        COMP_RHS(F, G, RHS, imax, jmax, delt, delx, dely);
        POISSON(P, RHS, imax, jmax, delx, dely, eps, itermax, omega, res);
        ADAP_UV(U, V, F, G, P, imax, jmax, delt, delx, dely);
        */

        COMP_DELT(delt, imax, jmax, delx, dely, U, V, Re, tau, problem);
        SETBCONDnew(U, V, P, FLAG, imax, jmax, wW, wE, wN, wS, problem);
        COMP_FG(U, V, F, G, FLAG, imax, jmax, delt, delx, dely, gx, gy, cgamma, Re);
        COMP_RHS(F, G, RHS, FLAG, imax, jmax, delt, delx, dely);
        POISSON(P, RHS, FLAG, imax, jmax, delx, dely, eps, itermax, omega, res, t);
        ADAP_UV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely);

        n += 1;
        if (n % 100 == 0) {
            write_data("finaldata/liddrivencavity200.txt", U, V, P, n);
        }
    }
    int t2 = time(0);
    int dur = t2 - t1;
    cout << dur << endl;

    write_data("finaldata/liddrivencavity200.txt", U, V, P, n);

    cout << "finished in " << n << " steps\n";

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




// in funktionen const flag