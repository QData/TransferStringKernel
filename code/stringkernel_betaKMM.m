%Implementation of Transfer String Kernel for TFBS prediction
%Original (k,m)-mismatch kernel code by Pavel P. Kuksa (http://filer2.org/)
%Modified code by Ritambhara Singh (rs3zz@virginia.edu)

function []=stringkernel_betaKMM(fastafile,kip,mip,ntrain)
%rs3zz:This code reads fasta file and returns kernel and label files
%rs3zz: Input format :fastafile=FASTA format input file, kip= k parameter value, mip= m parameter value,
%ntrain= number of training samples

seq = (fastaread( fastafile ));
nseq = length( seq );
fprintf(1,'# sequences: %d\n', nseq)

% Prepare sequence data.

% Map DNA sequences to strings over [0,3] alphabet
map('ACGT') = 0:3;
for i = 1:nseq
  % note that this removes all non-AGCT characters
  
  %rs3zz: added "upper" to handle lowercase characters inorder to fix the
  %k(i,i)=0 problem
  
  S{i} = map( regexprep(upper(seq(i).Sequence), '[^ACGT]','') );
end

% Now compute mismatch kernel.

k =kip; % k-mer length
m = mip; % maximum number of mismatches (m < k)
na = 4; % alphabet size

fprintf(1,'Computing (%d,%d)-mismatch kernel\n', k, m)
K = mismatchKernel( S, k, m, na ); 

K = double(K);
fprintf(1, 'done.\n')

fprintf(1,'Obtaining data labels...');
for i = 1:nseq
    str_labels{i} = regexprep(seq(i).Header,'^\S+\|','');
end
str_class_names = unique( str_labels );
nclass = length( str_class_names );
labels = -1*ones(nseq,1);

nsamples = ntrain;  % rs3zz: number of train samples
ntestsamples=nseq-ntrain; %rs3zz: number of test samples

%rs3zz:normalizing kernels:
Knorm=normalizeKernelMatrix(K);
Ktrainnorm=Knorm(1:nsamples,1:nsamples);
Ktestnorm=Knorm(nsamples+1:nseq,1:nsamples);

%rs3zz:change made to include varying training datasets

Kbetanorm=Knorm(1:nsamples,nsamples+1:nseq);
%rs3zz: printing training kernel matrix for SVM
kernel=fopen(strcat(fastafile,'.',num2str(k),'.',num2str(m),'.',num2str(ntrain),'.TRAINKERNEL.txt'),'a');
for p=1:nsamples
        if strcmp(str_labels{p},'-1')==1
            class_ind = -1;
        else
            class_ind = 1;
        end
        labels(p)=class_ind;
        fprintf(kernel,'%d ',labels(p));
        for col1=1:nsamples
            fprintf(kernel,'%d:%f ',col1,Ktrainnorm(p,col1));
        end
        %rs3zz: change in code for libsvm format
        %fprintf(kernel,'\n');
        %rs3zz: change in code for svm_light format
        fprintf(kernel,'# %d\n',p);
end
fclose(kernel);

%rs3zz: printing testing kernel matrix for SVM
kernel=fopen(strcat(fastafile,'.',num2str(k),'.',num2str(m),'.',num2str(ntrain),'.TESTKERNEL.txt'),'a');
lab=fopen(strcat(fastafile,'.LABELS.txt'),'w');
row=1;
for p=nsamples+1:nseq
        if strcmp(str_labels{p},'-1')==1
            class_ind = -1;
        else
            class_ind = 1;
        end
        labels(row)=class_ind;
        fprintf(lab,'%d\n',labels(row));
        fprintf(kernel,'%d ',labels(row));
        for col1=1:nsamples
            fprintf(kernel,'%d:%f ',col1,Ktestnorm(row,col1));
        end
	%rs3zz: change in code for libsvm format
        %fprintf(kernel,'\n');
        %rs3zz: change in code for svm_light format
        fprintf(kernel,'#%d \n',row);
row=row+1;
end
fclose(kernel);
fclose(lab);

 
% %betaKMM starts

% 'calculating H=K...'    
H = Ktrainnorm;
H=(H+H')/2; %make the matrix symmetric (it isn't symmetric before because of bad precision)
% 'calculating f=kappa...' 
R3 = Kbetanorm;
f=(R3*ones(ntestsamples, 1));
f=-nsamples/ntestsamples*f;

regression=0;
% % % subject to...
% % % abs(sum(beta_i) - m) <= m*eps
% % % which is equivalent to A*beta <= b where A=[1,...1;-1,...,-1] and b=[m*(eps+1);m*(eps-1)]
eps = (sqrt(nsamples)-1)/sqrt(nsamples);
%eps=1000/sqrt(nsamples);
A=ones(1,nsamples);
A(2,:)=-ones(1,nsamples);
b=[nsamples*(eps+1); nsamples*(eps-1)];

Aeq = [];
beq = [];
% 0 <= beta_i <= 1000 for all i
LB = zeros(nsamples,1);
UB = ones(nsamples,1).*1000;

% X=QUADPROG(H,f,A,b,Aeq,beq,LB,UB) attempts to solve the quadratic programming problem:
%              min 0.5*x'*H*x + f'*x   
% subject to:  A*x <= b 
%              Aeq*x = beq
%              LB <= x <= UB 

% 'solving quadprog for betas...'
[beta,FVAL,EXITFLAG] = quadprog(H,f,A,b,Aeq,beq,LB,UB);
EXITFLAG;
if ((EXITFLAG==0 ))
    %&& (doingRealTraining==1))
    [beta,FVAL,EXITFLAG] = quadprog(H,f,A,b,Aeq,beq,LB,UB,beta,optimset('MaxIter',1e6));
    EXITFLAG;
end

if (regression==0)
% guarantee that all beta greater than 0
    threshold=0.01*abs(median(beta));
    beta(find(beta<threshold)) = threshold;
    sprintf('number of beta < %f: %d (0 is good)', threshold, length(find(beta<threshold)))
end

% End of beta calculation

%rs3zz: printing weighted training kernel matrix for SVM
kernel=fopen(strcat(fastafile,'.',num2str(k),'.',num2str(m),'.',num2str(ntrain),'.WEIGHTTRAINKERNEL.txt'),'a');
for p = 1:nsamples
    if strcmp( str_labels{p},'-1')==1
        class_ind = -1;
    else
        class_ind = 1;
    end
    labels(p)=class_ind;
    fprintf(kernel,'%d ',labels(p));
    %rs3zz: change in code for libsvm format
    fprintf(kernel,'cost:%f ',beta(p));
    for col3=1:nsamples
        fprintf(kernel,'%d:%f ',col3,Ktrainnorm(p,col3));
    end
    %rs3zz: change in code for libsvm format
    %fprintf(kernel,'\n');
    %rs3zz: change in code for svm_light format
    fprintf(kernel,'# %d\n',p);
end
fclose(kernel);

 %rs3zz: change in code for libsvm format 
%storing beta values separately to
%input as weight vecotr in libsvm 
% betas=fopen(strcat(fastafile,'.BETA.txt'),'a');
% for p = 1:nsamples
%     fprintf(betas,'%d\n',beta(p));
% end
% fclose(betas);

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % EXITFLAG:
% %       1  QUADPROG converged with a solution X.
% %       3  Change in objective function value smaller than the specified tolerance.
% %       4  Local minimizer found.
% %       0  Maximum number of iterations exceeded.
% %      -2  No feasible point found.
% %      -3  Problem is unbounded.
% %      -4  Current search direction is not a direction of descent; no further 
% %           progress can be made.
% %      -7  Magnitude of search direction became too small; no further progress can
% %           be made.

%exit;


