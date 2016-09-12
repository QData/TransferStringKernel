function K= mismatchKernel(S,k,m,na)
%%MISMATCHKERNEL(S,k,m,na) computes mismatch(k,m) kernel for strings S
%  over alphabet of size na
%  O(n k^{m+1}) algorithm (independent of alphabet size)
%  K=mismatchKernel(S,k,m,na) computes mismatch(k,m) kernel for
%    strings S (cell array) over alphabet [0:na-1].
%
%  INPUTS:
%    1. S - cell array of strings, each cell S{i} is
%    a double array (string) with elements in the range [0, na-1],
%    where na is the alphabet size
%    2. k - k-mer length (eg. k=5)
%    3. m - max. number of mismatches (eg. m=1)     
%    3. na - alphabet size (eg. na=20)
%  OUTPUTS:  
%    1. K - kernel matrix (double matrix of size |S|-by-|S|)
%
%  Example:
%
%    S={[1 0 2 3 2 1],[1 1 2 0 1 1]}; % some sequences over [0:3] alphabet
%    K=mismatchKernel(S,5,1,4);       % mismatch(5,1) kernel
%

num_strings=length(S);

% Precompute mismatch weights.
w=zeros(min(2*m,k),1);
w(1)=1;
for i=1:m
 w(1)=w(1)+nchoosek(k,i)*(na-1)^i;
end

if (m==1)
 w(2) = na;
 w(3) = 2;
end

if (m==2)
 w(2) = 1 + k*(na-1) + (k-1)*(na-1)^2;
 w(3) = 1 + 2*(k-1)*(na-1) + (na-1)^2;
 w(4) = 6*(na-1);
 w(5) = 6;
end
if (m==3)
 w(2)=1+k*(na-1)+nchoosek(k,2)*(na-1)^2+nchoosek(k-1,2)*(na-1)^3;
 w(3)=1+k*(na-1)+(nchoosek(k,2)+nchoosek(k-2,2))*(na-1)^2+(k-2)*(na-1)^3;
 w(4)=1+3*(na-1)+(6*(k-3)+3)*(na-1)^2+(na-1)^3;
 w(5)=2+nchoosek(k-1,2)*(na-1)+12*(na-1)^2;
 w(6)=9+10*(na-1)+(na-1)^2;
 w(7)=nchoosek(6,3);
end


fprintf(1,'w=[ %s]\n', sprintf('%1.0f ',w));

% Extract k-mers.
[sxy resgroup nf]=mexExtractKmer(S,k);
nf=double(nf);
% Compute mismatch kernel.

for i=0:k-1
 comb{i+1}.comb=nchoosek([1:k],k-i);
end

nchoosekmat=zeros(k,k);

for i=k:-1:1
  for j=1:i
    nchoosekmat(i,j)=nchoosek(i,j);
  end
end

K=zeros(num_strings,num_strings);

for i=0:min(2*m,k)
 Ks(i+1).K=zeros(num_strings,num_strings);
end

for i=0:min(2*m,k-1)
  comb1=comb{i+1}.comb;
  num_iter=size(comb1,1);
  i
  for j=1:num_iter
    sxy1=sxy(:,comb1(j,:));
    
    % Compute spectrum kernel.
    regroup=mexcntsrtna(sxy1, k-i, na);
    sxy1=sxy1(regroup+1,:);
    resgroup1=resgroup(regroup+1);
    Ks(i+1).K=Ks(i+1).K+mexextractuniquena(sxy1, resgroup1, k-i, num_strings,na);
  end
end

num_max_mismatches=min(2*m,k);

for i=1:min(num_max_mismatches,k-1)%k-1
  for j=0:i-1
    Ks(i+1).K=Ks(i+1).K-nchoosekmat(k-j,i-j)*Ks(j+1).K;
  end
end
if (num_max_mismatches==k)
  Ks(k+1).K=nf*nf';
  for i=1:k
    Ks(k+1).K=Ks(k+1).K-Ks(i).K;
  end
end

%K=zeros(num_strings,num_strings,'int32');
for i=0:min(2*m,k)
  K=K+w(i+1)*Ks(i+1).K;
end
clear Ks
