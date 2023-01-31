n=100
a=speye(n,n)+sprandn(n,n,0.5/sqrt(n));
%A=full(a)
%[L,U]=lu(A)
[l,u]=ilu(a ,struct('type','ilutp','droptol',1e-2) );
figure(1),spy(a)
figure(2),spy(l,'+'),hold on,spy(u,'x'),hold off
iluErr=norm( a-(l*u) ,Inf)
fillFac=(nnz(l)+nnz(u))/nnz(a)
condesta=condest(a)
condestula=condest(u\(l\a))

