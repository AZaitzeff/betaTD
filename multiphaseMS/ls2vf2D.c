#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

double H(double x);
double IH(double x);
double IIH(double x);
int wrap(int j, int n);

/*************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{

  /* **********************************************************************
   * MATLAB Interface:
   * **********************************************************************
   * vf = ls2vf2D(x,y,lev,Z)
   * **********************************************************************
   * This routine returns a volume-of-fluid style representation for an
   * interface from an input of a level-set style representation.
   * **********************************************************************
   * INPUT:
   *  x = x-coords of pixels at which level-set representation is given.
   *      One dimensional array.
   *  y = y-coords of pixels at which level-set representation is given.
   *      One dimensional array.
   *  lev = Level set values at the x & y coordinates. 1D array.
   *  Z = Workspace variable; should be all -1's. Size should be the same
   *      as the 2D computational grid.
   * **********************************************************************
   */

  int *x, *y;
  int i, j, k, ind, xk, yk, E, W, N, S, Npix, m, n, test;
  double *lev, *Z, *vf;
  double mx, my, val;
  
  x = (int *) mxGetData(prhs[0]); /* Array of x coords. */
  y = (int *) mxGetData(prhs[1]); /* Array of y coords. */
  lev = (double *) mxGetData(prhs[2]); /* Level set values. */
  Z = (double *) mxGetData(prhs[3]); /* Workspace. */

  Npix = mxGetM(prhs[0]);  /* Number of pixels currently in the subset. */
  n = (int) mxGetN(prhs[3]); /* Dimension of grid. Assumed square. */
    
  /* Allocate memory and get pointer for output variables. */
  plhs[0] = mxCreateDoubleMatrix(Npix,1,mxREAL);	
  vf = mxGetPr(plhs[0]);
  
  /* Populate Z with level set values: */
  for (k=0;k<Npix;k++) { /* Loop over pixels. */
    xk = x[k] - 1;
    yk = y[k] - 1;
    ind = n*xk + yk;
    Z[ind] = lev[k];
    /* printf("%f\n",lev[k]); */
  }
  
  /* Real action: */
  for (k=0;k<Npix;k++){ /* Loop over pixels. */
    xk = x[k] - 1;
    yk = y[k] - 1;
    val = lev[k];
    ind = n*xk + yk;

    
    E = n*wrap(xk+1,n) + yk;
    W = n*wrap(xk-1,n) + yk;
    N = n*xk + wrap(yk+1,n);
    S = n*xk + wrap(yk-1,n);
        
    mx = 0.5*(Z[E]-Z[W]);
    my = 0.5*(Z[N]-Z[S]);
    
    test = 0;
    if (fabs(mx)<1e-6) test = 1;
    if (fabs(my)<1e-6) {
      test = 2;
      if (fabs(mx)<1e-6) test = 3;
    }
    
    if (test==0) {
      vf[k] = (IIH(mx+my+val) - IIH(-mx+my+val) - IIH(mx-my+val) + IIH(-mx-my+val)) / (mx*my);
    }
    
    if (test==1) {
        vf[k] = 2*( IH(my+val) - IH(val-my) ) / my;
    }
    
    if (test==2) {
        vf[k] = 2*( IH(mx+val) - IH(val-mx) ) / mx;
    }
    
    if (test==3) {
        vf[k] = 4*H(val);
    }
    
    vf[k] = vf[k] / 4.0;
    
    if ( max(max(max(max(Z[E],Z[W]),Z[N]),Z[S]),val) < 0 ) vf[k] = 0;
    if ( min(min(min(min(Z[E],Z[W]),Z[N]),Z[S]),val) > 0 ) vf[k] = 1.0;
    
  } /* for k */

  /* Clean up Z: */
  for (k=0;k<Npix;k++){
    xk = x[k]-1;
    yk = y[k]-1;
    ind = n*xk + yk;
    Z[ind] = -1.0;
  }

}
double H(double x) 
{
  if (x <= 0) return 0.0;
  else return 1.0;
}
  
double IH(double x) 
{
  if (x<=0) return 0.0;
  else return x;
}
  
double IIH(double x) 
{
  if (x<=0) return 0.0;
  else return 0.5*x*x;
}

int wrap(int j, int n)
{
  if (j<0) return (n-1);
  else {
    if (j>n-1) return 0;
  }
  return j;
}
