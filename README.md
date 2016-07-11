# Transfer String Kernel for Cross-Context Transcription Factor Binding Site (TFBS) Prediction 

Running the Matlab Code:

1. Download code folder and open in Matlab
2. Compile all the C files

>>mex mexcntsrtna.c
>>mex mexEtractKmer.c
>>mex mexextractuniqena.c
>>mex mexextractuniquena_int.c

3. Run TSK code

>>stringkernel_betaKMM('example.fasta',10,3,500)

(Format: stringkernel_betaKMM(fastafile,k parameter ,m parametr ,number of training samples)

4. This will generate the following kernel files:
