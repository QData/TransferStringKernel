### Transfer String Kernel for Cross-Context DNA-Protein Binding Prediction 

Running the Matlab Code:

1. Download code folder and open in Matlab
2. Compile all the C files
```matlab
mex mexcntsrtna.c
mex mexEtractKmer.c
mex mexextractuniqena.c
mex mexextractuniquena_int.c
```

Run TSK code
```matlab
stringkernel_betaKMM('example.fasta',10,3,1000)
```
(stringkernel_betaKMM(FASTA format file, k parameter , m parameter , number of training samples)

>Note: example.fasta consists of 2000 sequences, with following format:
```
>chr1:465738-478736 (position)|-1 (label)   [header]
TACGATCGATCGATCGATCGATCA.....ATCGCTCGAT (100bp length DNA sequences)
```

This will generate the following kernel files for input into SVM package SVMLight (commented code modifications available for LIBSVM format output in the code):

```
example.fasta.10.3.1000.TESTKERNEL.txt : Test Kernel file
example.fasta.10.3.1000.TRAINKERNEL.txt : Train Kernel file (Simple SK)
example.fasta.10.3.1000.WEIGHTTRAINKERNEL.txt : Train Kernel file with weights (TSK)
example.fasta.LABELS.txt : File containing true labels for the testing
```

***

#### Running SVM Classifier

We use the [SVMLight package](http://svmlight.joachims.org/) for implementing SVM classifier for TSK :

Once the kernel files are generated, use following commands to train/test the svm
>Note: These are just example commands, hyperparameter tuning of C parameter (-c) maybe be required to choose the best value.

Training:
```
svm_learn -c 1 example.fasta.10.3.1000.WEIGHTTRAINKERNEL.txt model.tmp
```
Classification:
```
svm_classify -f 1 example.fasta.10.3.1000.TESTKERNEL.txt model.tmp example.fasta.10.3.1000.WEIGHTPRED.txt
```

>Note: stringKernel_betaKMM.m also contains (commented) code to print out kernels in format for [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) package.
