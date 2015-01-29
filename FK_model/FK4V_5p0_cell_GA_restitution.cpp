#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <iostream>

#define NDIM 21 //number of parameters, used to be 23
#define DT 0.01
#define S1limit 10
#define BCL 500.

#define stim_dur 3.0
#define stim_amp 3

double heav(double u);

using namespace std;

double heav(double u) {
  if (u > 0.0)
    return 1.0;
  else
    return 0.0;
}

int main(int argc, char *argv[]) {

  // Simulation parameters
  long int n;
  double t, stim_time;
  int N_stim = round(stim_dur / DT);
  int N_BCL = round(BCL / DT);
  int numS1 = 0; 
  double TOTALTIME = S1limit * BCL + 1000.0;

  // State variables
  double du, dv, dw, dd, tso;
  double u = 0.0, uP = 0.0, v = 1.0, w = 1.0, d = 0.0, tempu;
  double xfi, xso, xsi, xstim;

  // Parameter vector
  vector<double> para(NDIM);
  
  //----------------KNOW PARAMETERS------------------//
  double uo, um, una;
  //---------------UNKNOWN PARAMETERS-------------------//
  double uc, uv, uw, ud, tvm, tvp, twm, twp, tsp, tsm, ucsi; 
  double xk, td, to, tsoa, tsob, uso, xtso, tsi, D, tvmm;

  ifstream input_file("xfinal.dat", ios::in);
  ofstream output_APD("FK4V_5p0_cell_GA_restitution_danParams.dat", ios::out);

  // Read parameters from file
  if(!input_file || !output_APD) {
    cerr << "open error!!" << endl;
    exit(1);
  }
  for(int i = 0; i < NDIM; i++) {
    input_file >> para[i];
  }
  input_file.close();
  
  // Restitution protocol parameters
  int DI[] = {200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 50, 40, 30, 25, 20, 15, 10, 5};
  // double DI[] = {254.0690, 222.4222, 182.5915, 152.8019, 123.0218, 103.2425, 83.4569, 43.7630, 24.1447, 14.3405};
  int N_DI;
  int size = sizeof(DI)/sizeof(DI[0]);

  
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
  for(int i = 0; i < size; i++) {

    t = 0.0;
    n = 0;
    numS1 = 0;
    N_DI = DI[i] / DT;
    stim_time = 0;

    // APD parameters
    // double u_max = -1.0, u_low = 0.04966, 
    double u_max = -1.0, u_low = 0.0;
    double u_90, obj_APD;
    double APD, APD_start = 0.0;

    while (t < TOTALTIME) {

      // Stimlulus Current
      if (n == stim_time) {
        numS1++;
        // cout << "STIM: n = " << n << " stim_time = " << stim_time << endl;
      }
      if (n - stim_time <= N_stim && n - stim_time > 0)
        xstim = stim_amp;
      else 
        xstim = 0.0;

      // Model Functions
      tso = tsoa + (tsob - tsoa) * (1 + tanh((u - uso) * xtso))/2.0;

      dv = (1 - heav(u - una)) * (1 - v) / ((1 - heav(u - uv)) * tvm + tvmm * heav(u - uv)) - heav(u - una) * v / tvp;
      dw = (1 - heav(u - uw)) * (1 - w) / twm - heav(u - uw) * w / twp;
      dd = ((1 - heav(u - ud)) / tsm + heav(u - ud) / tsp) * ((1 + tanh(xk * (u - ucsi))) / 2.0 - d);

      v = v + DT * dv;  //fast gate
      w = w + DT * dw;  //slow gate
      d = d + DT * dd;

      // Currents
      xfi = -v * heav(u - una) * (u - una) * (um - u) / td;
      xso = (u - uo) * (1 - heav(u - uc)) / to + heav(u - uc) / tso;
      xsi = -w * d / tsi;

      // du_dt
      du = - (xfi + xso + xsi - xstim);
      u = u + DT * du;

      // Determine time of next stimulus using APD_90
      if (n == stim_time - 1) {
        u_max = -1.0;
      }
      if (u > u_max) {
        u_max = u;
        u_90 = u_max - (u_max - u_low) * 0.90;
        // u_90 = 0.1;
        APD_start = t;
      }
      if (u <= u_90 && uP > u_90 && n - stim_time > 50) {
        APD = t - APD_start;
        cout << "time = " << t << " APD_start = " << APD_start << " APD = " << APD << " u_90 = " << u_90 << endl;
        if (numS1 < S1limit)
          stim_time = (numS1 + 1) * N_BCL;
        else if (numS1 == S1limit) {
          stim_time = n + N_DI; 
          cout << "last stim at: " << stim_time << endl;
        }
      }

      uP = u;
      
      // Output
      if (n % 10 == 0) {
        // output_voltage << u << " ";
      }

      n++; 
      t = n * DT;

    } // End of time loop
    
    obj_APD = 19.7238132356 * log(DI[i]) + 237.7581328626; //LA S2 rest
    output_APD << DI[i] + APD << "\t" << DI[i] << "\t" << APD << "\t" << obj_APD << endl;
    cout << "DI = " << DI[i] << " APD = " << APD << " obj_APD = " << obj_APD << endl;
    cout << "----------------------------------------" << endl;
    
  } // End s2 loop

}
