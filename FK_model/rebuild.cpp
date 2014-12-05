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
#define TOTALTIME 2400
#define DT 0.01
#define DX 0.02
#define SPACING_INDEX 10
#define S1limit 5

#define xstimdur 3.0
#define xstimamp 30.0

double heav(double u);

using namespace std;

double heav(double u) {
  if (u > 0.0)
    return 1.0;
  else
    return 0.0;
}

int main() {

  double savedt = DT;
  int numt = TOTALTIME/DT;
  int istimdur = round(xstimdur/DT);
  int savepoints = TOTALTIME/savedt;
  int i, j, k;
  
  int ni = 0;
  int ns = 0;
  double preV1 = 0, preV2 = 0;

  int numS1 = 0; //count how many S1 stim was applied (have ~5 inbetween each S2 stim)
  float S1time = 400./DT; //CL(cycle length) between S1 stim
  int lastStimTime = 0;
  int lastPeakTime = 0;

  double dv, dw, dd, tso;
  vector<double> u(NX), v(NX), w(NX), d(NX), tempu(NX + 1, 0.0);
  vector<double> xfi(NX), xso(NX), xsi(NX), xstim(NX);
  vector<double> para(NDIM);
  vector<double> excitimes(NEXCITE), nexcs(NEXCITE);
  vector<double> v70(savepoints + 1), v30(savepoints + 1);

  // for(i = 0; i < NEXCITE; i++) {
  //   excitimes[i] = TOTALTIME;
  // }

  // for (i = 0; i < NEXCITE; i++) {
  //   nexcs[i] = int(excitimes[i] / DT);
  // }  
  
  //----------------KNOW PARAMETERS------------------//
  double uo, um, una;
  //---------------UNKNOWN PARAMETERS-------------------//
  double uc, uv, uw, ud, tvm, tvp, twm, twp, tsp, tsm, ucsi; 
  double xk, td, to, tsoa, tsob, uso, xtso, tsi, D, tvmm;

  ifstream input_file("xfinal.dat", ios::in);
  ofstream output_file("rebuild_voltage.dat", ios::out);

  if(!input_file || !output_file) {
    cerr << "open error!!" << endl;
    exit(1);
  }

  for(i = 0; i < NDIM; i++) {
    input_file >> para[i];
    // cout << para[i] << endl;
  }

  input_file.close();
  
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

  for(i = 0; i < numt; i++) {

    for(j = 0; j < NX; j++) {

      if (i >= S1time + lastStimTime && i <= S1time + lastStimTime + istimdur && j < 5 && numS1 < S1limit) {
        if (j==4 && i + 1 > S1time + istimdur) {
          numS1++;
          lastStimTime = i;
          printf("S1 @ %i\n", i);
        }
        xstim[j] = xstimamp;
      }
      else if (i >= nexcs[ni] && i <= nexcs[ni] + istimdur && j < 5) {
        if (numS1 >= S1limit) {
          if (j == 4 && i + 1 > nexcs[ni] + istimdur) {
            ni++;
            numS1 = 0;
            lastStimTime = i;
          }
          xstim[j] = xstimamp;   
        }
      }
      else {
        xstim[j] = 0.0;
      }

      // Stimlulus Current
      // if(i == stim_time){
      //   numS1++;
      //   printf("Stim at time %i\n", i);
      // }
      // if(i - stim_time <= istimdur && i - stim_time > 0 && j < 5)
      //   xstim[j] = xstimamp;
      // else 
      //   xstim[j] = xstimamp;

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

    for(k = 0; k < NX; k++) {

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

    if (i % 10 == 0) {
      for (k = 0; k < NX; k++) {
        output_file << u[k] << " ";
      }
      output_file << endl;
    }

  }
}
