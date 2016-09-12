function [ Knorm ] = normalizeKernelMatrix( K )
%rs3zz: Function to normalize the kernel matrix s.t the diagonal value=1
    Kdiag=diag(K);
    Kroot=sqrt(Kdiag);
    k=1./Kroot;
    kt=k';
    Knorm=K.*(k*kt);
end

