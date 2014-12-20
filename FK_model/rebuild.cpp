#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <iostream>

#define NDIM 21 //number of parameters, used to be 23
#define EPS 1e-18
#define NX 100 //number of divisions in cable ("cells") (use ~100 for 1D, 2 for 0D)
#define NEXCITE 10 //number of excitations from stimulus current (-1 for APPLY S1)
#define DT 0.01
#define DX 0.02
#define SPACING_INDEX 10
#define S1limit 5
#define BCL 400.

#define xstimdur 3.0
#define xstimamp 30.0
#define APD_rec_site 3

double heav(double u);

using namespace std;

double heav(double u) {
  if (u > 0.0)
    return 1.0;
  else
    return 0.0;
}

int main() {

  long int n;
  double t;
  int istimdur = round(xstimdur/DT);
  int i, j, k;
  double TOTALTIME = S1limit * BCL + 1000.0;
  
  double preV1 = 0, preV2 = 0;

  int numS1 = 0; //count how many S1 stim was applied (have ~5 inbetween each S2 stim)
  float S1time = BCL / DT; //CL(cycle length) between S1 stim
  double stim_time = BCL / DT;

  double dv, dw, dd, tso;
  vector<double> u(NX), uP(NX), v(NX), w(NX), d(NX), tempu(NX + 1, 0.0);
  vector<double> xfi(NX), xso(NX), xsi(NX), xstim(NX);
  vector<double> para(NDIM);
  
  //----------------KNOW PARAMETERS------------------//
  double uo, um, una;
  //---------------UNKNOWN PARAMETERS-------------------//
  double uc, uv, uw, ud, tvm, tvp, twm, twp, tsp, tsm, ucsi; 
  double xk, td, to, tsoa, tsob, uso, xtso, tsi, D, tvmm;

  ifstream input_file("xfinal.dat", ios::in);
  ofstream output_voltage("rebuild_voltage_myStim.dat", ios::out);
  ofstream output_APD("rebuild_APD_restitution_cell3.dat", ios::out);

  // Read parameters from file
  if(!input_file || !output_voltage || !output_APD) {
    cerr << "open error!!" << endl;
    exit(1);
  }
  for(i = 0; i < NDIM; i++) {
    input_file >> para[i];
  }
  input_file.close();

  // Restitution protocol parameters
  // int DI[] = {200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 50, 40, 30, 25, 20, 15, 10, 5};
  int DI[] = {254.0690, 222.4222, 182.5915, 152.8019, 123.0218, 103.2425, 83.4569, 43.7630, 24.1447, 14.3405};
  int N_DI;
  int size = sizeof(DI)/sizeof(DI[0]);

  for (j = 0; j < NX; j++) {
    u[j] = 0.0;
    v[j] = 1.0;
    w[j] = 1.0;
    d[j] = 0.0;
  }
  
  //tvmm=10;
  uo = 0.0;
  um = 1.0;
  una = 0.23;
  uc = para[0]; //threshold
  uv = para[1]; //fast gate threshold, determines whethere tvm or tvmm is active (chaos 8)
  uw = para[2]; //slow gate threshold
  ud = para[3]; //threshold
  tvm = para[4]; //controls minimum diastolic interval where CV occurs (chaos 8)
  tvp = para[5]; //fast gate closing time
  twm = para[6]; //slow gate opening time (changes APD shape?)
  twp = para[7]; //slow gate closing time (shifts APD up/down?)
  tsp = para[8]; //d-gate variables
  tsm = para[9];
  ucsi = para[10];
  xk = para[11]; //typically around 10
  td = para[12]; //fast current time variable, determines max CV
  to = para[13]; //ungated time constant
  tsoa = para[14]; //curve shape/APD, ungated time, adjusts DI
  tsob = para[15]; //ungated time. Easily adjusts DI, changes APD
  uso = para[16];
  xtso = para[17];
  tsi = para[18]; //slow current time variable, max APD
  D = para[19]; //related to density, mostly changes CV, but can effect everything
  tvmm = para[20]; //controls the steepness of the CV curve (chaos 8)

  // Restitution protocol loop
  for(i = 0; i < size; i++) {

    t = 0.0;
    n = 0;
    numS1 = 0;
    N_DI = DI[i] / DT;
    stim_time = S1time;
    // t_end = S1limit * 400.0 + 1000.0;

    // APD parameters
    double u_max = -1.0, u_low = 0.04966, u_90, obj_APD;
    double APD, APD_start;
    


    while (t < TOTALTIME) {

      for(j = 0; j < NX; j++) {

        // Stimlulus Current
        if (n == stim_time && j == 0) {
          numS1++;
          cout << "STIM: n = " << n << " stim_time = " << stim_time << endl;
        }
        if (n - stim_time <= istimdur && n - stim_time > 0 && j < 5)
          xstim[j] = xstimamp;
        else 
          xstim[j] = 0.0;

        // Model Functions
        tso = tsoa + (tsob - tsoa) * (1 + tanh((u[j] - uso) * xtso))/2.0;

        dv = (1 - heav(u[j] - una)) * (1 - v[j]) / ((1 - heav(u[j] - uv)) * tvm + tvmm * heav(u[j] - uv)) - heav(u[j] - una) * v[j] / tvp;
        dw = (1 - heav(u[j] - uw)) * (1 - w[j]) / twm - heav(u[j] - uw) * w[j] / twp;
        dd = ((1 - heav(u[j] - ud)) / tsm + heav(u[j] - ud) / tsp) * ((1 + tanh(xk * (u[j] - ucsi))) / 2.0 - d[j]);

        v[j] = v[j] + DT * dv;  //fast gate
        w[j] = w[j] + DT * dw;  //slow gate
        d[j] = d[j] + DT * dd;

        //----------currents---------//
        xfi[j] = -v[j] * heav(u[j] - una) * (u[j] - una) * (um - u[j]) / td;
        xso[j] = (u[j] - uo) * (1 - heav(u[j] - uc)) / to + heav(u[j] - uc) / tso;
        xsi[j] = -w[j] * d[j] / tsi;
      
      }

      for (k = 0; k < NX; k++) {
        tempu[k] = u[k];
      }

      // PDE Integration
      for (k = 0; k < NX; k++) {
        if (k != 0 && k != NX - 1) {
          u[k] = tempu[k] + DT * (D * (tempu[k+1] + tempu[k-1] - 2*tempu[k])/pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k] - xstim[k]));
        }
        else {
          if (k == 0) {
            preV2 = preV1;
            preV1 = u[k];
            u[k] = tempu[k] + DT * (D * (tempu[k + 1] - tempu[k]) / pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k] - xstim[k]));
          }
          else {
            u[k] = tempu[k] + DT * (D * (tempu[k - 1] - tempu[k]) / pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k] - xstim[k]));
          }
        }
      }

      //APD_90
      if (n == stim_time - 1) {
        u_max = -1.0;
        //actual_DI = t - APD_end;
      }
      if (u[APD_rec_site] > u_max) {
        u_max = u[APD_rec_site];
        u_90 = u_max - (u_max - u_low)*0.90;
        APD_start = t;
      }
      if (u[APD_rec_site] <= u_90 && uP[APD_rec_site] > u_90 && n - stim_time > 50) {
        cout << "time = " << t << " APD = " << APD << " u_90 = " << u_90 << endl;
        APD = t - APD_start;

        if (numS1 <= S1limit)
          stim_time = (numS1 + 1) * S1time;
        else if (numS1 == S1limit + 1) {
          stim_time = n + N_DI; 
          cout << "last stim at: " << stim_time << endl;
        }
        // cout << "NEXT STIM: " << stim_time << endl;
      }

      for (k = 0; k < NX; k++) {
        uP[k] = u[k];
      }

      // Output
      if (n % 10 == 0) {
        for (k = 0; k < NX; k++) {
          output_voltage << u[k] << " ";
        }
        output_voltage << endl;
      }

      n++; 
      t = n * DT;

    } // End of time loop
    
    obj_APD = 19.7238132356*log(DI[i]) + 237.7581328626; //LA S2 rest
    output_APD << DI[i] + APD << "\t" << DI[i] << "\t" << APD << "\t" << obj_APD << endl;
    cout << "DI = " << DI[i] << " APD = " << APD << " obj_APD = " << obj_APD << endl;
    cout << "--------------------------" << endl;
    
  } // End s2 loop

}
