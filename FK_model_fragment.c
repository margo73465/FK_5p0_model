	//Fenton Karma Model
	//dd,dv,and dw are temporary double precision variables
	//DT = 0.01, DX = 0.02	

  for (j = 0; j < NX; j++) {
    tso = tsoa + (tsob - tsoa) * (1 + tanh((u[j] - uso) * xtso)) / 2.0;
    
    dv = (1 - heav(u[j] - una)) * (1 - v[j]) / ((1 - heav(u[j] - uv)) * tvm + tvmm * heav(u[j] - uv)) - heav(u[j] - una) * v[j] / tvp;
    dw = (1 - heav(u[j] - uw)) * (1 - w[j])/twm - heav(u[j] - uw) * w[j] / twp;        
    dd = ((1 - heav(u[j] - ud)) / tsm + heav(u[j] - ud) / tsp) * ((1 + tanh(xk * (u[j] - ucsi))) / 2.0 - d[j]);

    v[j] = v[j] + DT * dv;  //fast gate
    w[j] = w[j] + DT * dw;  //slow gate
    d[j] = d[j] + DT * dd;
            
    //----------currents---------//
    xfi[j] = -v[j] * heav(u[j] - una) * (u[j] - una) * (um - u[j]) / td;
    xso[j] = (u[j] - uo) * (1 - heav(u[j] - uc)) / to + heav(u[j] - uc) / tso;
    xsi[j] = -w[j] * d[j] / tsi;
            
  }
   
	//Laplacian  
  for (k = 0; k < NX; k++) {
    if (k != 0 && k != NX-1) {
      u[k] = tempu[k] + DT * (D * (tempu[k + 1] + tempu[k - 1] - 2 * tempu[k]) / pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k]));
    }
    else {
      if (k == 0) {
        preV2 = preV1;
        preV1 = u[k];
        u[k] = tempu[k] + DT * (D * (tempu[k + 1] - tempu[k]) / pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k]));
      }
      else {
        u[k] = tempu[k] + DT * (D * (tempu[k - 1] - tempu[k]) / pow(DX,2.0) - (xfi[k] + xso[k] + xsi[k]));
      }
    }
            
