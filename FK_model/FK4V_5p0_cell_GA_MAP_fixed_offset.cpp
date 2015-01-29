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
  long int n = 0;
  double t = 0.0;
  int N_stim = round(stim_dur/DT);
  int N_BCL = round(BCL/DT);
  int numS1 = 0;
  double TOTALTIME = S1limit * BCL;
  
  // State variables
  double dv, dw, dd, du, tso;
  double u = 0.0, uP = 0.0, v = 1.0, w = 1.0, d = 0.0;
  double xfi, xso, xsi, xstim;
  
  // Parameter vector
  vector<double> para(NDIM);

  // MAP match-up parameters
  double u_max, max_time;

  // GA fitness function parameters
  double diff, obj_u, sum_diff = 0.0;
  
  //----------------KNOW PARAMETERS------------------//
  double uo, um, una;
  //---------------UNKNOWN PARAMETERS-------------------//
  double uc, uv, uw, ud, tvm, tvp, twm, twp, tsp, tsm, ucsi; 
  double xk, td, to, tsoa, tsob, uso, xtso, tsi, D, tvmm;

  // for(int i = 1; i < argc - 1; i++){
  //   para[i-1] = atof(argv[i]);
  //   cout << para[i-1] << endl;
  // }

  ifstream input_file("best_params_FK4V_5p0_cell_GA_MAP_fixed_offset.dat", ios::in);
  ifstream objective("scaled_MAP.dat", ios::in);
  // ofstream out_error(argv[argc - 1], ios::out);
  ofstream out_voltage("FK4V_5p0_cell_GA_MAP_fixed_offset_voltage.dat", ios::out);

  for (int i = 0; i < NDIM; i++) {
    input_file >> para[i];
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
  

  while (t < TOTALTIME) {

    // Stimlulus Current
    if(n % N_BCL == 0) {
      numS1++;
    }
    if(n % N_BCL <= N_stim) xstim = stim_amp;
    else xstim = 0.0;

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

    // if (numS1 == S1limit - 1 && t - BCL * numS1 > 3.1 && t - BCL * numS1 < 403.1) {
    if (numS1 == S1limit - 1 && n % N_BCL > 300 && n % N_BCL < 50000 && n % 100 == 0) {
      objective >> obj_u; 
      diff = (u - obj_u) * (u - obj_u);
      sum_diff += diff;
      out_voltage << t << "\t" << u << "\t" << obj_u << "\t" << endl;
    }

    uP = u;

    n++; 
    t = n * DT;

  } // End of time loop

  if(sum_diff != sum_diff) {
    sum_diff = 9999999999;
  }
  // } else if (sum_diff < 0.0001) {
  //   cout << "too small!" << endl;
  //   sum_diff = 0.0001;
  // }

  cout << sum_diff << endl;
  // out_error << sum_diff << endl;
}
