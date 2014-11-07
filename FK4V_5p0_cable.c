/******************************************************************************
Fenton Karma Model by Daniel Lombardo
dd,dv,and dw are temporary double precision variables
DT = 0.01, DX = 0.02	

NOTES ON REFACTORING: 
  dv = dv_dt, dw = dw_dt, dd = dd_dt
  tvp = tau_v_plus, tvm = tau_v1_minus, tvmm = tau_v2_minus
  twm = tau_w_minus, twp = tau_w_plus
  tsm = tau_s_minus, tsp = tau_s_plus, xk = k_s, ucsi = u_csi
  uv = u_v, uw = u_w, ud = u_d, una = u_na (appears in I_fi), uso = u_so, 
    um = u_m, uo = u_o, uc = u_c
  tso = tau_so, tsoa = tau_so1, tsob = tau_so2, xtso = k_so 
  td = tau_d, to = tau_o, tso = tau_so, tsi = tau_si
  xfi = I_fi, xso = I_so, xsi = I_si

  tau_w_minus = in my version this has a similar equation to tau_v_minus
  r_s = also exists in my version but not here

Version 5p0 - entirely based on code from Daniel, and incorporated into a
              set-up that resembles my other models 

******************************************************************************/
#include <math.h>
#include <stdio.h>

const double DT = 0.01;
const double DX = 0.02;
const double BCL = 400;
const double stim_dur = 3.0;
const double stim_amp = 30.0;
const int paceAP = 5;
const int NX = 100;
const int APD_rec_site = 30;

int heav (double x);

int main() {

  // Files
  FILE *output_V;
  FILE *output_CV;
  FILE *output_APD; 
  output_V = fopen("FK4V_5p0_danielParams_BCL400_danielStim_V.dat", "w");
  output_CV = fopen("FK4V_5p0_danielParams_BCL400_danielStim_CV.dat", "w");
  output_APD = fopen("FK4V_5p0_danielParams_BCL400_danielStim_APD.dat", "w");

  // Variables
  int m, p, q, r;
  double tau_v_minus, tau_so;
  double dw_dt, dv_dt, dd_dt, du_dt;
  double v[NX], w[NX], d[NX], u[NX], tempu[NX], V[NX], VP[NX];
  double I_fi[NX], I_so[NX], I_si[NX], I_stim, cell_stim;
  double t, t_end, stim_time;
  int i, j, k, N;
  long int n;

  int N_BCL = BCL / DT;
  int N_stim = stim_dur / DT;

  // Restitution protocol parameters
  int DI[] = {200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 50, 40, 30, 25, 20, 15, 10, 5};
  int N_DI;
  int size = sizeof(DI)/sizeof(DI[0]);
  int count = 0;

  // Restitution protocol loop
  for(i = 0; i < size; i++) {

    t = 0.0;
    n = 0;
    N = -1;
    N_DI = DI[i] / DT;
    stim_time = 0.0;
    t_end = paceAP * BCL + 1000.0;

    // APD & CV variables
    double V_max = -100.0, V_low = -80.26, V_90;
    double APD, APD_end, APD_start = 0.0, obj_APD;
    double CV_start = 0.0;

    // Initialize state variables
    for (j = 0; j < NX; j++) {
      v[j] = 1.0;
      w[j] = 1.0;
      d[j] = 0.0;
      u[j] = 0.0;
      tempu[j] = 0.0;
      V[j] = 0.0;
      VP[j] = 0.0;
    }

    // Reset parameters
    double u_c = 2.09451397;
    double u_v = 0.32561365;
    double u_w = 0.22588498;
    double u_d = 0.46559551;
    double tau_v1_minus = 60.41746603;
    double tau_v2_minus = 1296.24007618;
    double tau_v_plus = 1.86551512;
    double tau_w_minus = 107.74420226;
    double tau_w_plus = 381.28267887;
    double tau_s_plus = 1.41854242;
    double tau_s_minus = 0.32968256;
    double u_csi = 0.32841700;
    double k_s = 5.09453859;
    double tau_d = 0.05371307;
    double tau_o = 36.98063948;
    double tau_so1 = 22.13423089;
    double tau_so2 = 19.82898648;
    double u_so = 0.36248000;
    double k_so = 7.98870889;
    double tau_si = 41.03797687;
    double u_o = 0.0;
    double u_m = 1.0;
    double u_na = 0.23;
    double D = 0.00068595;

    // Time loop
    while (t < t_end) {

      // Stimlulus Current
      if(n == stim_time){
        N++;
        printf("Stim at time %f\n", t);
      }
      if(n - stim_time <= N_stim && n - stim_time > 0)
        I_stim = stim_amp;
      else 
        I_stim = 0.0;

      // Cable loop
      for (j = 0; j < NX; j++) {

        cell_stim = 0.0;
        if (j < 5)
          cell_stim = I_stim;

        m = heav(u[j] - u_v);
        p = heav(u[j] - u_na);
        q = heav(u[j] - u_w);
        r = heav(u[j] - u_d);

        tau_v_minus = ((1 - m) * tau_v1_minus + tau_v2_minus * m);
        tau_so = tau_so1 + (tau_so1 - tau_so2) * (1 + tanh((u[j] - u_so) * k_so)) / 2.0;
        
        // State Variables
        dv_dt = (1 - p) * (1 - v[j]) / tau_v_minus - p * v[j] / tau_v_plus;
        dw_dt = (1 - q) * (1 - w[j]) / tau_w_minus - q * w[j] / tau_w_plus;        
        dd_dt = ((1 - r) / tau_s_minus + r / tau_s_plus) * ((1 + tanh(k_s * (u[j] - u_csi))) / 2.0 - d[j]);

        v[j] = v[j] + DT * dv_dt;  //fast gate
        w[j] = w[j] + DT * dw_dt;  //slow gate
        d[j] = d[j] + DT * dd_dt;
                
        // Currents
        I_fi[j] = -v[j] * p * (u[j] - u_na) * (u_m - u[j]) / tau_d;
        I_so[j] = (u[j] - u_o) * (1 - heav(u[j] - u_c)) / tau_o + heav(u[j] - u_c) / tau_so;
        I_si[j] = -w[j] * d[j] / tau_si;

        du_dt = - (I_fi[j] + I_so[j] + I_si[j] - cell_stim);
        tempu[j] = tempu[j] + DT * du_dt;

      } // End cable loop
       
      // // Laplacian (PDE solver) 
      // for (k = 0; k < NX; k++) {
      //   if (k != 0 && k != NX - 1) {
      //     //u[k] = tempu[k] + DT * (D * (tempu[k + 1] + tempu[k - 1] - 2.0 * tempu[k]) / pow(DX,2.0) - (I_fi[j] + I_so[j] + I_si[j] - cell_stim));
      //     u[k] = tempu[k] + DT * (D * (tempu[k + 1] + tempu[k - 1] - 2.0 * tempu[k]) / pow(DX,2.0));
      //   }
      //   else {
      //     if (k == 0) {
      //       //preV2 = preV1;
      //       //preV1 = u[k];
      //       //u[k] = tempu[k] + DT * (D * (tempu[k + 1] - tempu[k]) / pow(DX,2.0) - (I_fi[j] + I_so[j] + I_si[j] - cell_stim));
      //       u[k] = tempu[k] + DT * (D * (tempu[k + 1] - tempu[k]) / pow(DX,2.0));
      //     }
      //     else {
      //       //u[k] = tempu[k] + DT * (D * (tempu[k - 1] - tempu[k]) / pow(DX,2.0) - (I_fi[j] + I_so[j] + I_si[j] - cell_stim));
      //       u[k] = tempu[k] + DT * (D * (tempu[k - 1] - tempu[k]) / pow(DX,2.0));
      //     }
      //   }
      //   V[k] = 85.7 * u[k] - 84.0;
      // } // End PDE loop

      // PDE all inner grid points
      int cc, y;
      for (cc = 0; cc < 5; cc++) {
        for (y = 1; y < NX-1; y++)  {
          u[y] = DT/(DX * DX * 5) * D * (tempu[y + 1] + tempu[y - 1] - 2.0 * tempu[y]) + tempu[y];
        }
      
        // boundaries
        u[0] = u[2];
        u[NX - 1] = u[NX - 3];
      
        // update voltage in each cell to new value from PDE
        for (y = 0; y < NX; y++) {
          tempu[y] = u[y];
          V[y] = 85.7*u[y] - 84.0;
          //dV[y] = V[y] - VP[y];
        }
      }


      // CV
      if (N == paceAP) {
        if (V[19]>=-40. && VP[19]<-40.) CV_start = t;
        if (V[79]>=-40. && VP[79]<-40.) {
          //output_CV = fopen(output_CV_filename, "w");
          printf("time = %f\t t-CV_start = %f\t DX = %f\n", t, t - CV_start, DX);
          fprintf(output_CV,"%.2f\t%d\t%f\n", t, DI[i], (60.0 * DX / (t - CV_start) * 1000.0));
          //fclose(output_CV);
        }
      }

      //APD_90
      if (n == stim_time - 1){
        //V_low = V[APD_rec_site];  
        V_max = -100.0;
        //actual_DI = t - APD_end;
        //printf("DI = %d, APD reset, t = %f\n", DI[i], t);
      }
      if (V[APD_rec_site] > V_max) {
        V_max = V[APD_rec_site];
        V_90 = V_max - (V_max - V_low)*0.90;
        APD_start = t;
        //printf("APD_start = %f \t V_max = %f \t V_90 = %f\n",APD_start,V_max,V_90);
      }
      if (V[APD_rec_site] <= V_90 && VP[APD_rec_site] > V_90 && n - stim_time > 50){
        if(N < paceAP)
          stim_time = (N + 1) * BCL / DT;
        else if(N == paceAP){
          stim_time = n + N_DI; 
        }
        APD = t - APD_start;
        APD_end = t;
        count++;
        //printf("APD = %f\t time = %f\n", APD, t);
      }

      for (k = 0; k < NX; k++) {
        VP[k] = V[k];
      }

      if (n % 100 == 0) {
        for (j = 0; j < NX; j++) {
          fprintf(output_V, "%f\t", u[j]);
        }
        fprintf(output_V, "\n");
      }

      n++;
      t = n * DT;

    } // End time loop

    printf("APD now = %f\n", APD);
    obj_APD = 19.7238132356*log(DI[i]) + 237.7581328626; //LA S2 rest
    fprintf(output_APD, "%f\t%d\t%f\t%f\n", DI[i] + APD, DI[i], APD, obj_APD);
    
  } // End restitution loop

  fclose(output_V);
  fclose(output_CV);
  fclose(output_APD);

  return 0;

} // End main         

int heav (double x) {
  if (x >= 0) 
    return 1;
  else
    return 0;
}
