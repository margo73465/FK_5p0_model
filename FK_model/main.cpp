//
//  main.cpp
//  conduction1D
//
//  Created by Haicen on 1/5/13.
//  Copyright (c) 2013 Haicen. All rights reserved.
//

//if nan in apd, probably need more time
//if map not matching what it used to, check multiplier in ratio and main input
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include "nr3.h"

#define NDIM 21 //number of parameters, used to be 23
#define EPS 1e-18
#define NX 100 //number of divisions in cable ("cells") (use ~100 for 1D, 2 for 0D)
#define NEXCITE 10 //number of excitations from stimulus current (-1 for APPLY S1)
#define TOTALTIME 25000
#define DT 0.01
#define DX 0.02
#define SPACING_INDEX 10
#define S1limit 5

#define xstimdur 3.0
#define xstimamp 30.0

Doub heav(Doub u);

using namespace std;
Doub heav(Doub u) //Heaviside function
{
  if (u>0.0)
    return 1.0;
  else
    return 0.0;
}

struct solver_wave
{
  Doub savedt;
  int numt;
  bool success;
  int savepoints;
  VecInt nexcs;
  VecDoub DI;
  VecDoub u,v,w,d,xfi,xso,xsi;
  Doub dv,dw,dd,tso;
  VecDoub para,excitimes;
  VecDoub v70;
  VecDoub v30;
  
  Int istimdur;
  VecDoub xstim;
  
  //----------------KNOW PARAMETERS------------------//
  Doub uo, um, una;
  //---------------UNKNOWN PARAMETERS-------------------//
  Doub uc,uv,uw,ud,tvm,tvp,twm,twp,tsp,tsm,ucsi,xk,td,to,tsoa,tsob,uso,xtso,tsi,D,tvmm;
  
  solver_wave(VecDoub ppara,Doub ssavedt,VecDoub dDI,VecDoub eexcitimes):para(ppara),savedt(ssavedt),DI(dDI),excitimes(eexcitimes),u(NX+1,0.0),v(NX+1,0.0),w(NX+1,0.0),d(NX+1,0.0),xfi(NX+1,0.0),xso(NX+1,0.0),xsi(NX+1,0.0),xstim(NX+1,0.0),nexcs(NEXCITE+1,0),v70(100),v30(100)
  {
    //tvmm=10;
    uo=0.0;
    um=1.0;
    una=0.23;
    //t=0.0;
    
    savepoints=int(TOTALTIME/savedt);
    for (int j=0; j<NX; j++) {
      u[j]=0;
      v[j]=1;
      w[j]=1;
      d[j]=0;
      
    }

    numt=int(TOTALTIME/DT);
    
    istimdur=round(xstimdur/DT);
    
    uc=para[0]; //threshold
    uv=para[1]; //fast gate threshold, determines whethere tvm or tvmm is active (chaos 8)
    uw=para[2]; //slow gate threshold
    ud=para[3]; //threshold
    tvm=para[4]; //controls minimum diastolic interval where CV occurs (chaos 8)
    tvp=para[5]; //fast gate closing time
    twm=para[6]; //slow gate opening time (changes APD shape?)
    twp=para[7]; //slow gate closing time (shifts APD up/down?)
    tsp=para[8]; //d-gate variables
    tsm=para[9];
    ucsi=para[10];
    xk=para[11]; //typically around 10
    td=para[12]; //fast current time variable, determines max CV
    to=para[13]; //ungated time constant
    tsoa=para[14]; //curve shape/APD, ungated time, adjusts DI
    tsob=para[15]; //ungated time. Easily adjusts DI, changes APD
    uso=para[16];
    xtso=para[17];
    tsi=para[18]; //slow current time variable, max APD
    D=para[19]; //related to density, mostly changes CV, but can effect everything
    tvmm=para[20]; //controls the steepness of the CV curve (chaos 8)
    
    for (int i=0; i<NEXCITE; i++) {
      nexcs[i]=int(excitimes[i]/DT);
    }
    
    v70.resize(savepoints+1);
    v30.resize(savepoints+1);

    success=true;
  }
  void solving();
};
void solver_wave::solving()
{
  double maxVoltage[NEXCITE+1];
  int i,j,k;
  int ni=0;
  int ns=0;
  double preV1 =0, preV2=0;
  VecDoub tempu(NX+1,0.0);
  
  int numS1 = 0; //count how many S1 stim was applied (have ~5 inbetween each S2 stim)
  float S1time = 400./DT; //CL(cycle length) between S1 stim
  int lastStimTime = 0;
  int lastPeakTime = 0;
 
  
  FILE *fv = fopen("testv.txt","w+");
  //FILE *fxfi = fopen("xxfi.txt","w+");
  //FILE *fxsi = fopen("xxsi.txt","w+");
  //FILE *fxso = fopen("xxso.txt","w+");
  //FILE *fw = fopen("w.txt","w+");
  //FILE *fd = fopen("d.txt","w+");
  
  for (int init=0;init<NEXCITE;init++){
    maxVoltage[init]=0.0;
  }
  
  for(i=0;i<numt;i++)
  {

    for(j=0;j<NX;j++)
    {
      if (i >= S1time + lastStimTime && i<= S1time+lastStimTime+istimdur && j<5 && numS1 < S1limit) {
        if (j==4 && i+1 > S1time+istimdur) {
          numS1++;
          lastStimTime = i;
          //printf("S1 @ %i\n", i);
        }
        
        xstim[j] = xstimamp;
        
      }
      else if(i >= nexcs[ni] && i<= nexcs[ni]+istimdur && j<5) {
        if (numS1 >= S1limit) {
          if(j==4 && i+1 > nexcs[ni] + istimdur) {
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
      
      tso=tsoa+(tsob-tsoa)*(1+tanh((u[j]-uso)*xtso))/2.0;
      dv=(1-heav(u[j]-una))*(1-v[j])/((1-heav(u[j]-uv))*tvm + tvmm*heav(u[j]-uv)) - heav(u[j]-una)*v[j]/tvp;
      
      dw=(1-heav(u[j]-uw))*(1-w[j])/twm - heav(u[j]-uw)*w[j]/twp;
      
      dd=((1-heav(u[j]-ud))/tsm + heav(u[j]-ud)/tsp)*((1+tanh(xk*(u[j]-ucsi)))/2.0-d[j]);
      
      
      //printf("%f, %f, %f \n", dv,dw,dd);
      
      
      v[j]=v[j]+DT*dv;  //fast gate
      w[j]=w[j]+DT*dw;  //slow gate
      d[j]=d[j]+DT*dd;
      
      //----------currents---------//
      xfi[j]=-v[j]*heav(u[j]-una)*(u[j]-una)*(um-u[j])/td;
      xso[j]=(u[j]-uo)*(1-heav(u[j]-uc))/to + heav(u[j]-uc)/tso;
      xsi[j]=-w[j]*d[j]/tsi;
      
      //printf("%f \n", xfi[j]+xso[j]+xsi[j]);
      
    }
    
    //check max voltage (used for finding excitimes)
    int tempCheck = 3;
    if(S1limit<2) { //APPLY S1
      tempCheck = 0;
    }
    if(numS1>tempCheck) { //APPLY S1
      if(preV1 > preV2 && u[0] < preV1) { //max slope
        if(preV1 > maxVoltage[ni]) {
          maxVoltage[ni] = preV1;
          lastPeakTime = i;
          //printf("max volt %f @ %i\n", maxVoltage[ni],i);
        }
      }
    }
    
    if (i==round(ns*savedt/DT)&&ns<savepoints)
    {
      v70[ns] = u[70];
      v30[ns] = u[30];
      /*
       fprintf(fv, "%f \n", u[70]);
       fprintf(fxfi, "%4.15f \n", xfi[70]);
       fprintf(fxso, "%f \n", xso[70]);
       fprintf(fxsi, "%f \n", xsi[70]);
       fprintf(fw, "%f \n", w[70]);
       fprintf(fd, "%f \n", tsm);
       */
      if(i%10 == 1) {
        fprintf(fv, "%f \n", u[70]);
      }
      
      if (ns>0) {
        if (ns>lastStimTime*DT/savedt+100 && ni<NEXCITE) {
          float baseline = 0.1*maxVoltage[ni];
          //printf("%f \n", baseline);
          if(tempu[0]>baseline && u[0]<=baseline && numS1 >= S1limit)
          {
            //printf("%i \n", i - lastStimTime);
            if ((tempu[0]-baseline)>(baseline-u[0]))
            {
              nexcs[ni] = int(ns*savedt/DT) + int((DI[ni])/DT); //APPLY S1
            }
            else
            {
              nexcs[ni] = int((ns-1)*savedt/DT) + int((DI[ni])/DT);
            }
          }
        }
      }
      
      ns=ns+1;
      
      //printf("%i \n",ns);
    }
    
    for (k=0; k<NX; k++) {
      tempu[k]=u[k];
    }
    
    for(k=0;k<NX;k++)
    {
      
      if (k!=0 && k!=NX-1) {
        u[k] = tempu[k] + DT * (D * (tempu[k+1] + tempu[k-1] - 2*tempu[k])/pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k] - xstim[k]));
      }
      else {
        if (k==0) {
          preV2 = preV1;
          preV1 = u[k];
          u[k]=tempu[k]+DT*(D*(tempu[k+1]-tempu[k])/pow(DX,2.0)-(xfi[k]+xso[k]+xsi[k]-xstim[k]));
        }
        else {
          u[k]=tempu[k]+DT*(D*(tempu[k-1]-tempu[k])/pow(DX,2.0)-(xfi[k]+xso[k]+xsi[k]-xstim[k]));
        }
      }
      
    }
  }
  
  fclose(fv);
}

struct Funk
{
  int length; //from main(), represents number of data points (length of file)
  Doub savedt; //from main(), savedt = 1
  VecDoub data,DI; //actual experimental data
  VecDoub excitimes; //time of all excitations
  VecDoub para;
  int savepoints; //number of saved data points
  Doub error,error2,error3; //have three different sources of error (AP, APD, CV)
  int maxindex; //number excitations recorded for a certain "cell" or block
  int maxindex2;
  VecDoub timedata; //exp measured times (0,1,2,3,4,etc...)
  VecDoub maxslope; //record time of maximum slope of AP
  VecDoub maxslope2;
  VecInt maxmarker,maxmarker2; //mark the last recorded time of each excitation
  VecInt endmax, endmax2; //marking the end of a curve (goes with maxmarker/maxmarker2)
  VecDoub x_DI; //record diastolic interval and action potential duration
  VecDoub y_APD;
  VecDoub repol; //repolarization times
  
  int siteNum,siteNum2; //for site 70 and 30 respectfully
  
  //Constructor
  Funk(VecDoub ddata,VecDoub ttimedata,int llength,Doub ssavedt,VecDoub dDI,VecDoub eexcitimes):data(ddata),timedata(ttimedata),length(llength),savedt(ssavedt),DI(dDI),excitimes(eexcitimes),maxslope(NEXCITE,0.0),maxslope2(NEXCITE,0.0),maxmarker(NEXCITE,0),maxmarker2(NEXCITE,0),endmax(NEXCITE,0),endmax2(NEXCITE,0),x_DI(NEXCITE,0.0),y_APD(NEXCITE,0.0),repol(NEXCITE,0.0){
    savepoints=int(TOTALTIME/savedt);
  }
  
  Doub operator()(VecDoub para) //input parameters for fit, solve the model, check the error
  {
    solver_wave sol2(para,savedt,DI,excitimes); //take parameters, get simulation results
    sol2.solving();
    
    siteNum = 70; //usually 70, make both zero for zero d case
    siteNum2 = 30; //usually 30
    
    for(int ii=0;ii<NEXCITE;ii++)
    {
      //initial all excitation variables
      maxmarker[ii]=0; //time of slope change (max AP)
      maxmarker2[ii]=0;
      endmax[ii]=0;
      endmax[ii]=0;
      repol[ii]=0;
      maxslope[ii]=0;
      maxslope2[ii]=0;
      x_DI[ii]=0; //have a different DI and APD for each of the 8 stimulus currents
      y_APD[ii]=0;
    }
    
    maxindex=0;
    maxindex2=0;
    error=0.0;
    error2=0.0;
    error3=0.0;
    
    int i=0;
    int j=0;
    
    VecDoub slope2(savepoints-1,0.0); //record slope of AP curves
    VecDoub slope3(savepoints-1,0.0);
    
    int numS1 = 0; //count S1 stimuli
    int numS12 = 0;
    //int S1limit = 1;
    int lastS1Time = 0;
    int lastS1Time2 = 0;
    //int S1repol[NEXCITE]; //repolarization time of S1 stimuli right before S2 stimuli
  
    
    FILE *fmshape=fopen("mcurveshape.txt","w+");//same as fshape, but formatted for mathematica
    for (i=0; i<round(savepoints); i+=10) {
      //fprintf(fshape,"%8.4f,%8.6f\n",i*savedt,sol2.voltage[i][70]); //output AP vs. Time of cell 70 to curveshape
      fprintf(fmshape,"%8.6f \n",sol2.v70[i]);
    }
    //fclose(fshape);
    fclose(fmshape);
    
    //all results make use of the data from block 30 and 70 of the 100 piece heart cable
    //before computing actual error, fist check some physical properties of the fit
    for (i=0; i<savepoints-1; i++) {
      
      slope2[i]=sol2.v30[i+1]-sol2.v30[i]; //if data alright, record slope of action potential
      slope3[i]=sol2.v70[i+1]-sol2.v70[i];
    }
    
    for (i=0; i<savepoints-2; i++) {
      if (slope3[i]>0.0 && slope3[i+1]<=-0.0 && sol2.v70[i+1]>0.1) { //if slope changes sign
        //if (slope3[i]>0.0 && slope3[i+1]<=-0.0) {
        if (maxindex > NEXCITE) { //max index should never be greater than number of excitations
          printf("WRONG2, too many excitations %i\n",maxindex);
          return(1000);
        }
        
        if (maxindex >= 1) { //if we have already recorded at least one excitation
          //APPLY S1
          if (numS1 == 0) { //looking for first S1 stim
            if (i + 1 - maxmarker[maxindex - 1] < 100 / savedt) {
              if (sol2.v70[i + 1] > sol2.v70[maxmarker[maxindex - 1]]) {
                maxmarker[maxindex - 1] = i + 1; //still from the last excitation, so don't increment index
                //just move up the time for maxmarker
              }
            }
            else {
              numS1++;
              lastS1Time = i + 1;
              endmax[maxindex-1] = i+1;
            }
          }
          else if(numS1 < S1limit) { //now onto the next few S1 stimuli
            //add check for still on same stim
            if(i+1 - lastS1Time < 100/savedt) {
              if(sol2.v70[i+1]>sol2.v70[lastS1Time]) {
                lastS1Time = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              numS1++;
              lastS1Time = i+1;
            }
          }
          else { //either stil on last S1 stim or on an S2 stim
            if(i + 1 - lastS1Time < 100/savedt) {  //check for still on same stim
              if(sol2.v70[i + 1] > sol2.v70[lastS1Time]) {
                lastS1Time = i + 1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              maxmarker[maxindex] = i+1; //maxindex and maxmarker should be zero in this case
              maxindex++; //now have one excitation, marked at time i+1
              numS1 = 0;
            }
          }
          
        }
        else //looking for first excitation
        {
          //APPLY S1
          if(numS1 == 0) {
            numS1++;
            lastS1Time = i+1;
          }
          else if(numS1 < S1limit) {
            
            if(i+1 - lastS1Time < 100/savedt) {
              if(sol2.v70[i+1]>sol2.v70[lastS1Time]) {
                lastS1Time = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              numS1++;
              lastS1Time = i+1;
            }
            
          }
          else { //add check for still on same stim
            if(i+1 - lastS1Time < 100/savedt) {
              if(sol2.v70[i+1]>sol2.v70[lastS1Time]) {
                lastS1Time = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              maxmarker[maxindex] = i+1; //maxindex and maxmarker should be zero in this case
              maxindex++; //now have one excitation, marked at time i+1
              numS1 = 0;
            }
          }
          
        }
      }
      //if ((slope2[i]>0.0)&&(slope2[i+1]<=-0.0)&&maxindex2<=NEXCITE) { //same as above with other "cell"
      if (slope2[i]>0.0 && slope2[i+1]<=-0.0 && sol2.v30[i+1]>0.1) {
        if (maxindex2>NEXCITE) {
          printf("WRONG3, too many excitations\n");
          return(1000);
        }
        if (maxindex2>=1) {
          if (numS12 == 0) { //looking for first S1 stim
            if (i+1-maxmarker2[maxindex2-1]<100/savedt) {
              if (sol2.v30[i+1]>sol2.v30[maxmarker2[maxindex2-1]]) {
                maxmarker2[maxindex2-1]=i+1; //still from the last excitation, so don't increment index
                //just move up the time for maxmarker
              }
            }
            else {
              numS12++;
              lastS1Time2 = i+1;
              endmax2[maxindex2-1]=i+1;
            }
          }
          else if(numS12<S1limit) { //now onto the next few S1 stimuli
            //add check for still on same stim
            if(i+1 - lastS1Time2 < 100/savedt) {
              if(sol2.v30[i+1]>sol2.v30[lastS1Time2]) {
                lastS1Time2 = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              numS12++;
              lastS1Time2 = i+1;
            }
          }
          else { //add check for still on same stim
            if(i+1 - lastS1Time2 < 100/savedt) {
              if(sol2.v30[i+1]>sol2.v30[lastS1Time2]) {
                lastS1Time2 = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              maxmarker2[maxindex2] = i+1; //maxindex and maxmarker should be zero in this case
              maxindex2++; //now have one excitation, marked at time i+1
              numS12 = 0;
            }
          }
          
        }
        else
        {
          if(numS12 == 0) {
            numS12++;
            lastS1Time2 = i+1;
          }
          else if(numS12 < S1limit) {
            
            if(i+1 - lastS1Time2 < 100/savedt) {
              if(sol2.v30[i+1]>sol2.v30[lastS1Time2]) {
                lastS1Time2 = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              numS12++;
              lastS1Time2 = i+1;
            }
            
          }
          else { //add check for still on same stim
            if(i+1 - lastS1Time2 < 100/savedt) {
              if(sol2.v30[i+1]>sol2.v30[lastS1Time2]) {
                lastS1Time2 = i+1; //same S1 stim, just move up lastStimTime
              }
            }
            else {
              maxmarker2[maxindex2] = i+1; //maxindex and maxmarker should be zero in this case
              maxindex2++; //now have one excitation, marked at time i+1
              numS12 = 0;
            }
          }
   
        }

      }
    }
    //  printf("maxindex=%d,maxindex2=%d\n",maxindex,maxindex2);
    if (maxindex!=NEXCITE||maxindex2!=NEXCITE) { //both 30 and 70 should have recorded NEXCITE excitations
      if(maxindex<NEXCITE||maxindex2<NEXCITE){
        printf("WRONG4\n");
        printf("maxindex1 = %i,maxindex2 = %i\n",maxindex,maxindex2);
        return(1000);
      }
    }
    
    //====================//
    //          MAP       //
    //====================//
    
    VecDoub curveshape(savepoints,0.0); //record the shape of the AP for each excitation
    int stopindex,stopdata; //stopindex marks the time that a particular AP curve ends
    VecDoub datashape(savepoints,0.0);
    VecDoub normtime(length,0.0);
    Doub ratio=0.0;
    int dataSpacing;
    FILE *mModelShape = fopen("mModelShape.txt","w+");
    FILE *mDataShape = fopen("mDataShape.txt","w+");
    
    for (j=0; j<NEXCITE-1; j++) {
      if (maxmarker[j+1]-maxmarker[j] <= 1) { //time of each excitation should occur in the correct order
        printf("WRONG5\n");
        return(1000);
      }
      
      stopindex=0;
      stopdata=0;
      for (i = maxmarker[j]; i < endmax[j]; i++) {
        curveshape[i - maxmarker[j]] = sol2.v70[i]; //curveshape is AP
      }
      for (i=0+10./DT; i<endmax[j]-maxmarker[j]-1; i++) {
        if (curveshape[i]>0.1&&curveshape[i+1]<=0.1) { //AP curve has ended
          stopindex=i; //mark the time
          //printf("%i \n", stopindex);
          break;
        }
      }
      for (i=0; i<length-1; i++) {
        if (data[i]>0.1&&data[i+1]<=0.1) {
          stopdata=i; //same idea as curveshape, but with exp data
          //printf("%i \n", i);
          break;
        }
      }
      //printf("stopindex=%d,stopdata=%d\n",stopindex,stopdata);
      float datat = 1.0; //doesn't really matter, but make sure its the same multiplier as in main (time data)
      
      ratio=(float(stopindex*savedt))/(stopdata*datat); //convert stopindex to units of data, find ratio
      datashape[0]=data[0];
      
      for (i=0; i<length; i++) {
        normtime[i]=timedata[i]*ratio; //normalize the time values so its easier to compare with curveshape
      }
      error=error+fabs(datashape[0]-curveshape[0])/(0.01*curveshape[0]); //(exp-model)/sigma
      
      //based on the length of the data file, choose ~30 or ~10 evenly placed points to actually consider the error
      dataSpacing = round(stopdata/SPACING_INDEX);
      int dataIndex = 0;
      
      int jstart=0;
      
      for (i = 1; i < stopindex + 1; i++) { //for all times of a single curve
        for (int k = jstart; k < length - 1; k++) { //go over the length of the data file
          
          if (normtime[k] < i * savedt && normtime[k + 1] >= i * savedt) {
            
            datashape[i]=(data[k+1]-data[k])/(normtime[k+1]-normtime[k])*(i*savedt-normtime[k])+data[k];
            
            dataIndex++;
            
            if(dataIndex > dataSpacing) {
              dataIndex = 0;
            }
            
            jstart=k+1;
            
            fprintf(mModelShape,"%8.6f \n",curveshape[i]);
            fprintf(mDataShape,"%8.6f \n",datashape[i]);
            
            break;
          }
        }
      }
     
    }
    
    fclose(mDataShape);
    fclose(mModelShape);
    
    //=====================//
    //    APD Analysis     //
    //=====================//
    
    Doub baseline=0.0; //define line between DI and APD
    double S1repol[NEXCITE];
    float b;
    
    for (j=0;j<NEXCITE;j++) //find all maximum
    {
      baseline=0.1*sol2.v70[maxmarker[j]];
      for (int k=maxmarker[j]-10./DT;k<endmax[j]; k++) { //start at maxmarker - 5 ms
        if (sol2.v70[k]>baseline && sol2.v70[k-1]<baseline){
          b = sol2.v70[k]-slope3[k-1]*k;
          maxslope[j]=(baseline-b)/slope3[k-1];
          break;
        }
      }
    }
    
    for (j=0;j<NEXCITE;j++)
    {
      baseline=0.1*sol2.v30[maxmarker2[j]];
      
      for (int k=maxmarker2[j]-10./DT; k < endmax2[j]; k++) {
        
        if (sol2.v30[k]>=baseline && sol2.v30[k-1]<baseline){
          b = sol2.v30[k] - slope2[k-1]*k;
          maxslope2[j]=(baseline - b)/slope2[k-1];
          break;
        }
      }
    }
    
    for (j=1;j<NEXCITE;j++)
    {
      baseline=0.1*sol2.v70[maxmarker[j-1]];
      for (int k=maxmarker[j-1]+ 1./DT;k<endmax[j-1]-1./DT; k++) { //looking for repolarization times
        if((sol2.v70[k]>baseline)&&(sol2.v70[k+1]<=baseline)) //crossing baseline
        {
          b = sol2.v70[k+1] - slope3[k]*(k+1);
          repol[j-1] = (baseline - b)/slope3[k];
          break;
        }
      }
    }
    
    //now for the last excitation that we can measure
    baseline=0.1*sol2.v70[maxmarker[NEXCITE-1]]; //for the last possible diastolic interval
    
    for (int k=maxmarker[NEXCITE-1] + 1./DT; k < endmax[NEXCITE-1]- 1./DT; k++) {
      if(sol2.v70[k] > baseline && sol2.v70[k+1]<=baseline)
      {
        b = sol2.v70[k+1] - slope3[k]*(k+1);
        repol[j-1] = (baseline - b)/slope3[k];
        break;
      }
    }
    
    //find time that last S1 stim ends
    baseline=0.1;
    for (int k=0;k<maxslope[0]-1; k++) { //looking for repolarization times
      if((sol2.v70[k]>baseline)&&(sol2.v70[k+1]<=baseline)) //crossing baseline
        //going from APD to DI
      {
        b = sol2.v70[k+1] - slope3[k]*(k+1);
        S1repol[0] = (baseline - b)/slope3[k];
      }
    }
    
    for (j=1;j<NEXCITE;j++) //very similar to above
    {
      baseline=0.1;
      for (int k=maxmarker[j-1];k<maxslope[j]-1; k++) {
        if((sol2.v70[k]>baseline)&&(sol2.v70[k+1]<=baseline)) //crossing baseline
          //going from APD to DI
        {
          b = sol2.v70[k+1] - slope3[k]*(k+1);
          S1repol[j] = (baseline - b)/slope3[k];
        }
      }
    }
    
    for (j=0;j<NEXCITE;j++) //define the APD and DI
    {
      x_DI[j]=(maxslope[j]-S1repol[j])*savedt;
      
      y_APD[j]=(repol[j]-maxslope[j])*savedt;
    }
    
    double temp=0.0;
    FILE *f_APD=fopen("APD.txt","w+");
    
    for(i=0; i<NEXCITE; i++){
      fprintf(f_APD,"%8.4f %8.4f \n",x_DI[i],y_APD[i]);
    }
    
    fclose(f_APD);
    
    
    //=====================//
    //    CV analysis      //
    //=====================//
    temp=0.0;
    i=0;
    j=0;
    VecDoub speed(NEXCITE,0.0);
    FILE *f_speed=fopen("speed.txt","w+");
    
    for (j=0; j<NEXCITE; j++) {
      speed[j]=(abs(siteNum-siteNum2)*DX)/(0.001*savedt*abs(maxslope[j]-maxslope2[j])); //speed = dx/dt
      
    }
    
    for (i=0; i<NEXCITE;i++){
      fprintf(f_speed,"%8.4f %8.4f \n",x_DI[i],speed[i]);
    }
    
    fclose(f_speed);
    
    return(0);
  }
};

int main(int argc, char *argv[])
{
  VecDoub para(NDIM);
  double runningtime;
  clock_t begin,finish;
  Doub savedt = DT;
  
  ifstream infile1("xfinal.dat",ios::in);
  
  if(!infile1)
  {
    cerr<<"open error!!"<<endl;
    exit(1);
  }
  for(int i=0;i<NDIM;i++)
  {
    infile1>>para[i];
    cout << para[i] << endl;
  }
  infile1.close();
  
  int kmax=2000;
  VecDoub time(kmax,0.0),excitimes(NEXCITE,0.0),data(kmax,0.0),DI(NEXCITE,0.0);
  int length;
  
  for(int i=0;i<NEXCITE;i++) {
    excitimes[i]=TOTALTIME;
  }
  
  //load the experimental data
  int i=0;
  double maxMAP = -100.;
  double minMAP = 100000.;
  ifstream infile3("MAP.txt",ios::in);
  if(!infile3)
  {
    cerr<<"open error!"<<endl;
    exit(1);
  }
  float datat = 1.0;
  while(infile3.eof()==0)
  {
    time[i]=i*datat;
    infile3>>data[i];
    
    if (maxMAP<data[i])
      maxMAP = data[i];
    if (minMAP>data[i])
      minMAP = data[i];
    
    i++;
  }
  length=i-1;
  infile3.close();
  
  //normalize the MAP data
  for (int j=0;j<length;j++) {
    data[j] = (data[j]-minMAP)/(maxMAP-minMAP);
  }
  
  i=0;
  int j=0;
  
  ifstream infile5("DI.txt",ios::in);
  if(!infile5)
  {
    cerr<<"no DI data!"<<endl;
    exit(1);
  }
  while(infile5.eof()==0)
  {
    infile5>>DI[j];
    j++;
  }
  infile5.close();
  
  begin=clock();
  //solver_wave sol2(para,savedt,DI,excitimes); //take parameters, get simulation results
  //sol2.solving();
  Funk funk(data,time,length,savedt,DI,excitimes);
  funk(para);
  finish=clock();
  runningtime=(double)(finish-begin)/CLOCKS_PER_SEC;
  printf("the running time is %f s after the %d itertation.\n",runningtime,i);
  
  return 0;
}
