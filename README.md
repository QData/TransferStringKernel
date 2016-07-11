#### Transfer String Kernel for Cross-Context Transcription Factor Binding Site (TFBS) Prediction 

Running the Matlab Code:

1. Download code folder and open in Matlab
2. Compile all the C files

```matlab
mex mexcntsrtna.c
mex mexEtractKmer.c
mex mexextractuniqena.c
mex mexextractuniquena_int.c
```

3. Run TSK code

```matlab
stringkernel_betaKMM('example.fasta',10,3,1000)
```
>example.fasta consists of 2000 sequences 

(Format: stringkernel_betaKMM(FASTA format file, k parameter , m parameter , number of training samples)

4. This will generate the following kernel files for input into SVM package SVMLight (commented code modifications available for LIBSVM format output in the code):

+ example.fasta.10.3.1000.TESTKERNEL.txt : Test Kernel file
+ example.fasta.10.3.1000.TRAINKERNEL.txt : Train Kernel file (Simple SK)
+ example.fasta.10.3.1000.WEIGHTTRAINKERNEL.txt : Train Kernel file with weights (TSK)
+ example.fasta.LABELS.txt : File containing true labels for the testing


