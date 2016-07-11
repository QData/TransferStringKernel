#include "math.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
     inputs:
        k-mer features list sx
        features labels g
        feature length k
        number of strings nStr
        alphabet size na
     outputs:
        kernel matrix
  */
  
  const mxArray *xData, *gData, *k_arg, *na_arg, *nstr_arg;
  int *sx1, *g;
  long int r, c, i, j, ucnt, startInd, endInd, j1, cu;
  int k, na, nStr;
  char same;
  int *outK, *curfeat, *ucnts, *updind;
  int old_g;

  xData = prhs[0]; /* features */
  gData = prhs[1]; /* feature labels */
  k_arg = prhs[2]; /* k */
  nstr_arg = prhs[3];
  na_arg = prhs[4];

  sx1 = (int*)mxGetPr(xData);
  g = (int*)mxGetPr(gData);
  r = mxGetM(xData);
  c = mxGetN(xData);
  k = (int)(mxGetScalar(k_arg));
  na = (int)mxGetScalar(na_arg);
  nStr = (int)(mxGetScalar(nstr_arg));
  curfeat = mxCalloc(k, sizeof(int));
  plhs[0] = mxCreateNumericMatrix(nStr, nStr, mxINT32_CLASS, mxREAL);
  outK = (int*)mxGetPr(plhs[0]);

  i = 0;
  ucnts = (int *)mxCalloc(nStr, sizeof(int));
  updind = (int *)mxCalloc(nStr, sizeof(int));
  ucnt = 0;
  while (i<r)
  {
    for (j=0;j<k;j++)
      curfeat[j]=sx1[i+j*r];
    same=1;
    for (j=0;j<k;j++)
      if (curfeat[j]!=sx1[i+j*r])
      {
        same=0;
        break;
      }
    startInd=i;
    while (same && i<r)
    {
      i++;
      if (i>=r) break;
      same=1;
      for (j = 0; j < k; ++j)
      {
        if (curfeat[j] != sx1[i+j*r])
        {
          same=0;
          break;
        }
      }
    }
    endInd= (i<r) ? (i-1):(r-1);
    ++ucnt;
    
    /*if ((long int)endInd-startInd+1>2*nStr)  
    {*/
    old_g = -10;
    cu = 0;
    for (j = startInd; j <= endInd; ++j)
    {
      if (old_g != g[j])
      {
         updind[cu++]=g[j];
         ucnts[g[j]]=0;
         old_g = g[j];
      }
      ucnts[g[j]]++;
    }
    /*
      for (j = 0; j < nStr; ++j) ucnts[j] = 0;
      for (j = startInd; j <= endInd; ++j) ucnts[g[j]]++;
      cu = 0;
      for ( j = 0; j < nStr; ++j)
        if (ucnts[j]>0)
        {
          updind[cu]=j; cu++;
        }
     */
     /*mexPrintf("cu=%d\n",cu);
     mexPrintf("n = %d\n",((long int)endInd-startInd+1));*/
      for ( j = 0; j < cu; ++j)
        for (j1 = 0; j1 < cu; ++j1)
          outK[updind[j]+updind[j1]*nStr]+=ucnts[updind[j]]*ucnts[updind[j1]];
    /*}
    else
    {
      mexPrintf("in here");
      for (j=startInd;j<=endInd;j++)
        for (j1=startInd;j1<=endInd;j1++)
          outK[ g[j]+nStr*g[j1] ]++;
    }*/
  }
  mxFree(ucnts); mxFree(updind); 
}
