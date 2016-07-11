#include <math.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* inputs:
         sx - list of k-mers
         k - k-mer length
         na - alphabet size
       outputs:
         sorting permutation
    */
    char isVerbose;
    int *sx1;
    int k, na;
    int r, c;
    int *sxc, *sxl, *bc, *cc, *regroup;
    int i, j;

    isVerbose = 0;
    sx1 = (int*)mxGetPr(prhs[0]);
    r = mxGetM(prhs[0]);
    c = mxGetN(prhs[0]);
    na = (int)mxGetScalar(prhs[2]);
    k = (int)(mxGetScalar(prhs[1]));
    if (isVerbose)
    {
        mexPrintf("na=%d\n",na);
        mexPrintf("r=%d, c=%d", r, c);
        mexPrintf("Computing spectrum kernel k=%d\n", k);
    }
    sxc = mxCalloc(na, sizeof(int));
    bc  = mxCalloc(na, sizeof(int));
    sxl = mxCalloc(r, sizeof(int));
    cc  = mxCalloc(r, sizeof(int));
    plhs[0] = mxCreateNumericMatrix(r, 1, mxINT32_CLASS, mxREAL);
    regroup=(int*)mxGetPr(plhs[0]);
    for (i = 0; i < r; ++i)
       regroup[i]=i;
    for (j = k-1; j >= 0; --j)
      {
            for ( i = 0; i < na; ++i)
               sxc[i] = 0;
            for (i = 0; i < r; ++i)
            {
                cc[i] = sx1[ regroup[i]+j*r ];
                sxc[ cc[i] ]++;
            }
            bc[0]=0;
            for (i = 1;i < na; ++i)
               bc[i] = bc[i-1] + sxc[i-1];
            for (i = 0; i < r; ++i)
               sxl[bc[ cc[i] ]++] = regroup[i];
            for (i = 0; i < r; ++i)
               regroup[i] = sxl[i];
      }
   mxFree(sxl); mxFree(cc); mxFree(sxc); mxFree(bc);

   return;
}
