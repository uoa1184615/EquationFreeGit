% test some parallel spmd aspects, AJR 12 Dec 2020
%echo on
n=2000
ord=4
spmd
%{ 
Create a spmd distribution over the processors in the second
dimension. Then create a distributed square array of random
numbers.
%}
thecodist=codistributor1d(2)
disp('A ='); A=randn(n,thecodist)

%{ 
The first alternative is what I have currently coded which
is to make a working array with extra dimensions first, then
the size of A, and create working arrays index by the first
dimension.  But the distribution is in the second dimension
of suma, not the third, so it seems that execution time here
is longer by having to re-organise the data in A across the
cpus.
%}
disp('suma ='); 
suma=zeros([ord size(A)],thecodist)
tic, suma(1,:,:)=A;
for p=2:ord, suma(p,:,:)=cumsum(suma(p-1,:,:),2); end
sumaTime=toc

%{ 
The second alternative is to make a working array with extra
dimensions last, and first the size of A, and create working
arrays index by the last dimension.  It seems that execution
time here is about halved because it should not have to
re-organise the data in A across the cpus.  However, when
n=5000 then the times are roughly the same??
%}
disp('asum ='); 
asum=zeros([size(A) ord],thecodist)
tic, asum(:,:,1)=A;
for p=2:ord, asum(:,:,p)=cumsum(asum(:,:,p-1),1); end
asumTime=toc

%{
check the difference between two codings is small.
%}
%differ=gather(squeeze(suma(ord,:,:)))-gather(asum(:,:,ord));
%normdiffer=max(abs(differ(:)))
end

%echo off
