#include "math.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int *sx1, *g;
        int k, nStr, na;
        int *curfeat, *ucnts, *updind;
        int cu;
        double *outK;
        char same;
        long int r, c, i, j, ucnt;
        long int startInd, endInd, j1;

        sx1 = (int *)mxGetPr(prhs[0]);
	g = (int *)mxGetPr(prhs[1]);
        k=(int)mxGetScalar(prhs[2]);
        nStr=(int)mxGetScalar(prhs[3]);
        na=(int)mxGetScalar(prhs[4]);
	r=mxGetM(prhs[0]);
        c=mxGetN(prhs[0]);
	/*mexPrintf("r=%d, c=%d", r, c);*/

	/*printf("Computing spectrum kernel k=%d\n", k);*/
        /*mexPrintf("Number of strings=%d\n", nStr);*/
	curfeat=mxCalloc(k, sizeof(int));
	plhs[0]=mxCreateNumericMatrix(nStr, nStr, mxDOUBLE_CLASS, mxREAL);
	outK=mxGetPr(plhs[0]);
	/*mexPrintf("Extracting unique k-mers...\n");*/
	i=0;
        ucnt=0;
	ucnts=(int *)mxCalloc(nStr, sizeof(int));
	updind=(int *)mxCalloc(nStr, sizeof(int));

	while (i<r)
	{
		for (j = 0; j < k; ++j)
			curfeat[j]=sx1[i+j*r];
		same=1;
		for (j = 0;j < k; ++j)
		{
			if (curfeat[j]!=sx1[i+j*r])
			{
				same=0;
				break;
			}
		}
		startInd=i;
		while (same && i<r)
		{
			i++;
			if (i>=r) break;
			same=1;
			for (j = 0; j < k; ++j)
			{
				if (curfeat[j]!=sx1[i+j*r])
				{
					same=0;
					break;
				}
			}
		}
		endInd= (i<r) ? (i-1):(r-1);
		ucnt++;
		/*mexPrintf("unique feature: %d:%d:%d\n", ucnt, startInd, endInd);*/
/*if (startInd-endInd+1>2*nStr)	
{*/
		for (j = 0; j < nStr; ++j)
                       ucnts[j]=0;
	        for (j = startInd; j <= endInd; ++j)
                       ucnts[g[j]]++;
                cu=0;
		for (j = 0; j < nStr; ++j)
		{
  		   if (ucnts[j]>0)
		   {
			updind[cu]=j;
                        ++cu;
		   }
		}
		for (j=0;j<cu;j++)
		  for (j1=0;j1<cu;j1++)
			outK[updind[j]+updind[j1]*nStr]+=(double)ucnts[updind[j]]*ucnts[updind[j1]];
/*}
else
{
		int ind1;
		for (j=startInd;j<=endInd;j++)
		{
			for (j1=startInd;j1<=endInd;j1++)
			{
			/*ind1=g[j]+nStr*g[j1];*/
			/*	outK[ g[j]+nStr*g[j1] ]++;
				
			}
		}

 }              */
		
		
	}
        mxFree(ucnts); mxFree(updind); 
	/*mexPrintf("Number of unique features: %d\n", ucnt);*/
}
