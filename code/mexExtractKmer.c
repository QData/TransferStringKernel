#include "math.h"
#include "mex.h"
#include <stdbool.h>



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*
	inputs:
		cell array of sequences
		parameter k
	outputs:
		extracted features
		feature labels
                number of features per sequence
	*/
        bool isVerbose=true;
	const mxArray *cellS = prhs[0];
	if (!mxIsCell(cellS))
	{
		mexErrMsgTxt("***Error: input is not a cell array!\n");
	}
	int nStr=mxGetNumberOfElements(cellS);	
	int k=(int)(mxGetScalar(prhs[1]));
        if (isVerbose)
        {
 	  mexPrintf("Number of strings:%d\n", nStr);
          mexPrintf("k=%d\n",k);	
        }
	mxArray *sPtr;double *sBuf; int sLen;
	int i, j, j1, c;
	int *sx, *g;
        plhs[2]=mxCreateNumericMatrix(nStr, 1, mxINT32_CLASS, mxREAL);
        int *nf=(int*)mxGetPr(plhs[2]);
	int sumLen=0, numF=0;
	for (i=0;i<nStr;i++)
	{
		sPtr=mxGetCell(cellS, i);
		sLen=mxGetNumberOfElements(sPtr);
		sumLen+=sLen;
                nf[i]=((sLen>=k)?(sLen-k+1):0);
		numF+=nf[i];
	}
        if (isVerbose)
        {
	  mexPrintf("numF=%d, sumLen=%d\n", numF, sumLen);
        }
	/* Here we assume that k < (mininum string length) */
	plhs[0]=mxCreateNumericMatrix(numF, k, mxINT32_CLASS, mxREAL);
	plhs[1]=mxCreateNumericMatrix(numF, 1, mxINT32_CLASS, mxREAL);
	sx=(int*)mxGetPr(plhs[0]); g=(int*)mxGetPr(plhs[1]);
	c=0;
	for (i=0;i<nStr;i++)
	{
		sPtr = mxGetCell(cellS, i);
		sLen = mxGetNumberOfElements(sPtr);
                sBuf = mxGetPr(sPtr);
		for (j=0;j<sLen-k+1;j++)
		{
			for (j1=0;j1<k;j1++)
			{
				sx[c+j1*numF]=(int)sBuf[j+j1];
			}
			g[c]=i;
			c++;
		}
	}
	mexPrintf("Number of features extracted=%d\n", c);
	
	return;
}
